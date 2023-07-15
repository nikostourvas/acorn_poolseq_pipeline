#!/bin/bash

mkdir -p ../results/frequencies

# declare variables
BAM=../results_bak/align_bak/
OUTPUT=../results/frequencies

grenedalf sync \
    --sam-path $BAM/*Q20.markdup.bam \
    --out-dir $OUTPUT \
    --file-prefix bamfiles_ \
    --log-file $OUTPUT/grenedalf_sync_log.txt \
    --allow-file-overwriting