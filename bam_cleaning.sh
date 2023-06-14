#!/bin/bash

# declare variables
IND=${1}
OUTPUT=../results/align/${IND}

# Re-sort the reads according to their genomic coordinates
samtools sort -@ 16 ${OUTPUT}.fixmate.sort.bam > ${OUTPUT}.fixmate.newsort.bam
# Quality filter 
samtools view -@ 16 -bh -q 20 ${OUTPUT}.fixmate.newsort.bam > ${OUTPUT}.sort.Q20.bam
# mark the duplicate reads as variant callers take this information into account
# here we don't use 'samtools rmdup' as it is no longer recommended to remove duplicates
# marking the duplicates is enough
samtools markdup -@ 16 ${OUTPUT}.sort.Q20.bam ${OUTPUT}.sort.Q20.markdup.bam
# index
samtools index -@ 16 ${OUTPUT}.sort.Q20.markdup.bam
# gather stats
samtools flagstat ${OUTPUT}.sort.Q20.markdup.bam > ${OUTPUT}.sort.Q20.markdup.flagstat
samtools idxstats ${OUTPUT}.sort.Q20.markdup.bam > ${OUTPUT}.sort.Q20.markdup.idxstat
samtools stats -@ 16 ${OUTPUT}.sort.Q20.markdup.bam > ${OUTPUT}.sort.Q20.markdup.stats
samtools depth -@ 16 -aa ${OUTPUT}.sort.Q20.markdup.bam > ${OUTPUT}.sort.Q20.markdup.depth

# remove unnecessary files
rm ${OUTPUT}.sort.bam ${OUTPUT}.sort.Q20.bam
