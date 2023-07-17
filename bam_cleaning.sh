#!/bin/bash

# declare variables
IND=${1}
OUTPUT=../results/align/${IND}

# Re-sort the reads according to their genomic coordinates NIKOS
samtools sort -@ 2 ${OUTPUT}.fixmate.sort.bam > ${OUTPUT}.fixmate.newsort.bam
rm ${OUTPUT}*.sort.bam
# Quality filter 
samtools view -@ 2 -bh -q 20 ${OUTPUT}.fixmate.newsort.bam > ${OUTPUT}.sort.Q20.bam
# mark (and remove) the duplicate reads as variant callers take this information into account
# here we don't use 'samtools rmdup' as it is no longer maintained. 
# 'markdup -r' should achieve the same as 'rmdup'
samtools markdup -r -s -@ 16 ${OUTPUT}.sort.Q20.bam ${OUTPUT}.sort.Q20.markdup.bam

# index
samtools index -@ 2 ${OUTPUT}.sort.Q20.markdup.bam
# gather stats
samtools flagstat ${OUTPUT}.sort.Q20.markdup.bam > ${OUTPUT}.sort.Q20.markdup.flagstat
samtools idxstats ${OUTPUT}.sort.Q20.markdup.bam > ${OUTPUT}.sort.Q20.markdup.idxstat
samtools stats -@ 2 ${OUTPUT}.sort.Q20.markdup.bam > ${OUTPUT}.sort.Q20.markdup.stats
    # the following is extremely IO intensive - find alternative
#samtools depth -@ 16 -aa ${OUTPUT}.sort.Q20.markdup.bam > ${OUTPUT}.sort.Q20.markdup.depth

# remove unnecessary files
rm ${OUTPUT}.sort.Q20.bam
