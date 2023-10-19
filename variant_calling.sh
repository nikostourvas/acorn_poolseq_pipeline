#!/bin/bash

# declare variables
BAM_LIST=/data/genetics_tmp/variant_calling_tmp_storage_all_pools/AllPoolBams_TechnicalDupesRm_FastStorage.txt
OUTDIR=/data/genetics_tmp/VCF_AllPools
REF=/mnt/reference/Qrob_PM1N.fa
CHUNK=$1
CHUNK_SHORT=$(basename ${CHUNK/.bed/})
THREADS=1

mkdir -p ${OUTDIR}

#Create a .txt file that contains the file names (without directory information or suffixes) found in the BAM_LIST.
#First, remove any possible older versions of this file

rm /data/genetics_tmp/variant_calling_tmp_storage_all_pools/SampleNaming_VCF_${CHUNK_SHORT}.txt &&

while read line; do 
SAMPLE_NAME=$(basename ${line} | cut -d "." -f 1);
printf "%s\n" "${SAMPLE_NAME}" >> /data/genetics_tmp/variant_calling_tmp_storage_all_pools/SampleNaming_VCF_${CHUNK_SHORT}.txt;
done < ${BAM_LIST}

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

samtools mpileup -B -q 20 -l ${CHUNK} -f ${REF} -b ${BAM_LIST} -o ${OUTDIR}/${CHUNK_SHORT}_samtools.mpileup \
	2> ${OUTDIR}/${CHUNK_SHORT}.mpileup.err &&

java -jar /usr/share/java/varscan.jar mpileup2cns ${OUTDIR}/${CHUNK_SHORT}_samtools.mpileup\
	--vcf-sample-list /data/genetics_tmp/variant_calling_tmp_storage_all_pools/SampleNaming_VCF_${CHUNK_SHORT}.txt \
	--min-coverage 30 \
	--min-var-freq 0.025 \
        --min-reads2 1 \
        --min-freq-for-hom 0.75 \
        --p-value 0.1 \
        --output-vcf 1 \
        2> ${OUTDIR}/${CHUNK_SHORT}.varScan.snpindel.err \
        	| bgzip --compress-level -1 2> ${OUTDIR}/${CHUNK_SHORT}_Gzipping.err \
		> ${OUTDIR}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz &&

# index vcfs
bcftools index ${OUTDIR}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz \
    --threads ${THREADS} \
    2> ${OUTDIR}/${CHUNK_SHORT}.bcftools_index.snpindel.vcf.err

# extract snps and save them in a separate compressed VCF
bcftools view -v snps --threads ${THREADS} \
    ${OUTDIR}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz \
    -Oz -o ${OUTDIR}/${CHUNK_SHORT}.varScan.snp.vcf.gz \
    2> ${OUTDIR}/${CHUNK_SHORT}.bcftools_view.snp.vcf.err \

# extract indels and save them in a separate compressed VCF
bcftools view -v indels --threads ${THREADS} \
    ${OUTDIR}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz \
    -Oz -o ${OUTDIR}/${CHUNK_SHORT}.varScan.indel.vcf.gz \
    2> ${OUTDIR}/${CHUNK_SHORT}.bcftools.indel.vcf.err

# remove redundant files for storage efficiency
rm ${OUTDIR}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz \
   ${OUTDIR}/${CHUNK_SHORT}.varScan.snpindel.vcf.gz.csi

# some scaffolds will have no snps/indels
# delete the empty vcf files originating from these scaffolds
# find ${OUTDIR}/${CHUNK_SHORT}*.gz -maxdepth 1 -type f -empty -print -delete
#Defunct, but useful for any necessary benchmarking in the future. 

# index newly created VCFs
bcftools index ${OUTDIR}/${CHUNK_SHORT}.varScan.snp.vcf.gz \
    --threads ${THREADS} \   
    2> ${OUTDIR}/${CHUNK_SHORT}.bcftools_index.snp.vcf.err

bcftools index ${OUTDIR}/${CHUNK_SHORT}.varScan.indel.vcf.gz \
    --threads ${THREADS} \
    2> ${OUTDIR}/${CHUNK_SHORT}.bcftools_index.indel.vcf.err &&

echo -e "All log files for ${CHUNK_SHORT/.bed/}\n\n#####\n\nSamtools mpileup\n\n">${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.mpileup.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nvarScan.snpindel\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools index\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.snpindel.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools view snp\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_view.snp.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools view indels\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools.indel.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools index snp\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.snp.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools index snp\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.indel.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log &&

rm ${OUTDIR}/${CHUNK_SHORT}*.err 
rm ${OUTDIR}/${CHUNK_SHORT}_samtools.mpileup
