#!/bin/bash

# Check FASTQ file integrity
# Generate a hash for each file
mkdir -p ../results
md5sum /mnt/data/untrimmed_fastq/*.fastq.gz >  ../results/hashList.txt &&

# make list of samples to allow for parallel analysis of multiple samples
for i in /mnt/AcornSeqdata/AdapterClipped_Batch2/Sample*/*_R1_clipped.fastq.gz
do
	echo $(basename ${i%_R*})
done > /mnt/AcornSeqdata/inds_Batch2.txt

# Fastp for read trimming and deduplication
parallel --verbose -j 42 \
	'bash fastp_dedup.sh  {}' :::: ../Acorn_SeqData/AdapterClipped_Batch2/fastqnames.txt 

# Create a unified report with MultiQC
multiqc ../results -o ../results

# Map reads to reference genome
# You can change number of cores per job by accessing mapping.sh script.
/usr/local/bin/bwa-mem2/bwa-mem2 index /mnt/reference/Qrob_PM1N_with_cp_mt.fa

parallel --verbose -j 90 \
	'bash map.sh {}' :::: ../AcornSeqdata/samplenames.txt

# Filter out contigs <500 kb

parallel --verbose -j 90 \
	'bash bam_filtering_no_markdup.sh {}' :::: /data/genetics_tmp/results/fastp_dedup_trim/samplenames_bySize.txt

# make list of all chromosomes & scaffolds for next step
bash Simplified_Chunks.sh ../reference/Qrob_PM1N.fa.fai 2000000

# Variant calling
parallel --verbose -j 80 \
	'bash /mnt/acorn_poolseq_pipeline/variant_calling.sh {}' :::: /data/genetics_tmp/REFERENCE/ChunkFiles/Locations_Of_Chunk_Beds.txt

# Merge the multiple small VCF files that were produced in the previous step
# into one large VCF
bash Concatenate_VCFs.sh ../results/[VCF_DIRECTORY] [PREFIX_CONCATENATED_VCF] [THREADS]

# Filter SNPs detected close to INDELs, as suggested by VarScan manual
bash snp_indel_rm.sh

# Filter out multi-allelic SNPs
bash vcf_biallelic.sh

# Create a unified final report with MultiQC
multiqc --force ../results -o ../results
