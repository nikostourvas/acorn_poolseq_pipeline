#Lars Littmann
#2023.08.23
#Call variants (SNPs and INDELs) from BAMs that contain INDIVIDUAL sequencing data (not pool sequenced).
#For the ACORN project, these individuals are 2x20 Quercus robur trees belonging to populations 305 and 312
#This script is designed to run in parallel for regions of the genome.
#In the case of the ACORN project, these regions are arbitrarily made 'Chunks' of roughly equal size.
#Much of this script was adapted from Nikolaos Tourvas' script 'variant_calling.sh', part of the acorn poolseq pipeline

BAM_IN=../results/align_Batch1
OUTDIR=../results/indVCF
REF=../reference/Qrob_PM1N.fa
CHUNK=$1
CHUNK_SHORT=$(basename ${CHUNK/.bed/})
THREADS=1

# Create a mpileup file for each genomic region and call snps & indels together
# Input: (i) Filtered BAM files, (ii) indexed reference genome
# Output: (i) compressed genomic region VCF files for SNPs,
#         (ii) compressed genomic region VCF files for INDELs

# Multiple instances of this script can be run concurrently with GNU parallel.
# The mpileup | varscan | bgzip pipeline requires 2 CPU threads to run
# efficiently. For this reason we assign 2 CPU threads to the downstream
# bcftools steps.

# samtools mpileup arguments:
# -B: disable Base Alignment Quality (BAQ) adjustment as recommended by VarScan
# publication (Kobolt et al., 2013)
# -q: minimum mapping quality
# -r: region where variant calling is performed
# -f: reference genome file

# varscan arguments:
# --vcf-sample-list: a list of sample names in the same order as BAM files,
# one per line
# --min-coverage: Minimum read depth at a position to make a call
# --min-var-freq: Minimum variant allele frequency threshold
# --min-reads2: Minimum supporting reads at a position to call variants
# --min-freq-for-hom: Minimum frequency to call homozygote
# --p-value: P-value threshold for variant calling. It is recommended to use
# less stringent values for INDELs than SNPs. Here we call INDELs & SNPs
# together so we use p=0.1 for both. We will filter SNPs at a later
# stage though with the script "snp_indel_rm.sh".
# --output-vcf: Set to 1, to produce VCF file instead of table of alleles

samtools mpileup -B -q 20 -l ${CHUNK} -f ${REF} ${BAM_IN}/*Pl1-?????.markdup.Q20.bam \
    2> ${OUTDIR}/${CHUNK_SHORT}.mpileup.err \
    | java -jar /usr/share/java/varscan.jar mpileup2cns \
            --vcf-sample-list inds \
            --min-coverage 5 \
            --min-var-freq 0.025 \
            --min-reads2 1 \
            --min-freq-for-hom 0.85 \
            --p-value 0.1 \
            --output-vcf 1 \
            2> ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.snpindel.err \
            | bgzip --compress-level -1 \
                                > ${OUTDIR}/${CHUNK}.Ind.varScan.snpindel.vcf.gz

# index vcfs
bcftools index ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.snpindel.vcf.gz \
    --threads ${THREADS} \
    2> ${OUTDIR}/${CHUNK_SHORT}.Ind.bcftools_index.snpindel.vcf.err

# extract snps and save them in a separate compressed VCF
bcftools view -v snps --threads ${THREADS} \
    ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.snpindel.vcf.gz \
    -Oz -o ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.snp.vcf.gz \
    2> ${OUTDIR}/${CHUNK_SHORT}.Ind.bcftools_view.snp.vcf.err \

# extract indels and save them in a separate compressed VCF
bcftools view -v indels --threads ${THREADS} \
    ${OUTDIR}/${CHUNK_SHORT}.IND.varScan.snpindel.vcf.gz \
    -Oz -o ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.indel.vcf.gz \
    2> ${OUTDIR}/${CHUNK_SHORT}.Ind.bcftools.indel.vcf.err

# remove redundant files for storage efficiency
rm ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.snpindel.vcf.gz \
   ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.snpindel.vcf.gz.csi

# some scaffolds will have no snps/indels
# delete the empty vcf files originating from these scaffolds
find ${OUTDIR}/${CHUNK_SHORT}*.gz -maxdepth 1 -type f -empty -print -delete

# index newly created VCFs
bcftools index ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.snp.vcf.gz \
    --threads ${THREADS} \
    2> ${OUTDIR}/${CHUNK_SHORT}.Ind.bcftools_index.snp.vcf.err

bcftools index ${OUTDIR}/${CHUNK_SHORT}.Ind.varScan.indel.vcf.gz \
    --threads ${THREADS} \
    2> ${OUTDIR}/${CHUNK_SHORT}.Ind.bcftools_index.indel.vcf.err
