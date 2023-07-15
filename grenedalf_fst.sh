#!/bin/bash

mkdir -p ../results/frequencies

# declare variables
OUTPUT=../results/frequencies
POOLSIZES=/mnt/data/pool_sizes

grenedalf fst \
    --window-type sliding \
    --window-sliding-width 1000 \
    --method unbiased-nei \
    --pool-sizes $POOLSIZES \
    --filter-sample-min-count 8 \
    --filter-sample-max-coverage 500 \
    --filter-sample-min-coverage 50 \
    --filter-total-only-biallelic-snps \
    --sync-path $OUTPUT/bamfiles_sync.sync \
    --out-dir $OUTPUT \
    --log-file $OUTPUT/grenedalf_fst_log.txt \
    --allow-file-overwriting \
    --threads 4