#!/bin/bash

mkdir -p ../results/frequencies

# declare variables
OUTPUT=../results/frequencies
SYNC=/mnt/workshop_backup/bamfiles_sync.sync.gz
POOLSIZES=/mnt/data/pool_sizes

grenedalf diversity \
    --window-type sliding \
    --window-sliding-width 1000 \
    --pool-sizes $POOLSIZES \
    --filter-sample-min-count 8 \
    --filter-sample-max-coverage 500 \
    --filter-sample-min-coverage 50 \
    --filter-total-only-biallelic-snps \
    --sync-path $SYNC \
    --out-dir $OUTPUT \
    --log-file $OUTPUT/grenedalf_div_log.txt \
    --allow-file-overwriting \
    --threads 4