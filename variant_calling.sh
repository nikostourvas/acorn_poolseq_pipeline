#!/bin/bash

# create output directory
mkdir -p ../results/VCF

# declare variables
BAM=../results/Test_VariantCall_Input
RESULTS=../results/VCF
REF=../reference/Qrob_PM1N.fa
CHUNK=$1
CHUNK_SHORT=$(basename ${CHUNK/.bed/})
THREADS=2

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

# bgzip is used inside this pipeline to avoid writing large VCFs to disk
# --compress-level: Default value -1 (minimal compression). This setting is
# very fast. This can be set from -1 up to 9. However, higher compression
# levels will cost in runtime.

# Note from testing with subsampled files: If not even a single read is mapped 
# to a scaffold then the resulting mpileup file will be empty and varscan will 
# wait forever for input. However this is quite unlikely with real data.

samtools mpileup -B -q 20 -l ${CHUNK} -f ${REF} ${BAM}/*Pl1-???.markdup.Q20.bam \
	2> ${RESULTS}/${CHUNK_SHORT}.mpileup.err \
    | java -jar /usr/share/java/varscan.jar mpileup2cns \
            --vcf-sample-list inds \
            --min-coverage 30 \
            --min-var-freq 0.025 \
            --min-reads2 1 \
            --min-freq-for-hom 0.75 \
            --p-value 0.1 \
            --output-vcf 1 \
            2> ${RESULTS}/${CHUNK_SHORT}.varScan.snpindel.err \
            | bgzip --compress-level -1 \
				> ${RESULTS}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz &&

# index vcfs
bcftools index ${RESULTS}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz \
    --threads ${THREADS} \
    2> ${RESULTS}/${CHUNK_SHORT}.bcftools_index.snpindel.vcf.err

# extract snps and save them in a separate compressed VCF
bcftools view -v snps --threads ${THREADS} \
    ${RESULTS}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz \
    -Oz -o ${RESULTS}/${CHUNK_SHORT}.varScan.snp.vcf.gz \
    2> ${RESULTS}/${CHUNK_SHORT}.bcftools_view.snp.vcf.err \

# extract indels and save them in a separate compressed VCF
bcftools view -v indels --threads ${THREADS} \
    ${RESULTS}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz \
    -Oz -o ${RESULTS}/${CHUNK_SHORT}.varScan.indel.vcf.gz \
    2> ${RESULTS}/${CHUNK_SHORT}.bcftools.indel.vcf.err

# remove redundant files for storage efficiency
# rm ${RESULTS}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz \
#   ${RESULTS}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz.csi

# some scaffolds will have no snps/indels
# delete the empty vcf files originating from these scaffolds
find ${RESULTS}/${CHUNk_SHORT}*.gz -maxdepth 1 -type f -empty -print -delete

# index newly created VCFs
bcftools index ${RESULTS}/${CHUNK_SHORT}.varScan.snp.vcf.gz \
    --threads ${THREADS} \   
    2> ${RESULTS}/${CHUNK_SHORT}.bcftools_index.snp.vcf.err

bcftools index ${RESULTS}/${CHUNK_SHORT}.varScan.indel.vcf.gz \
    --threads ${THREADS} \
    2> ${RESULTS}/${CHUNK_SHORT}.bcftools_index.indel.vcf.err
