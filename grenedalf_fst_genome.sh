#!/bin/bash

mkdir -p ../results/frequencies

# declare variables
OUTPUT=../results/frequencies
SYNC=../results/frequencies/simulate.sync.gz
POOLSIZES=pool_sizes_sim

grenedalf fst \
    --window-type genome \
    --method unbiased-nei \
    --pool-sizes $OUTPUT/$POOLSIZES \
    --filter-sample-min-count 8 \
    --filter-sample-max-coverage 500 \
    --filter-sample-min-coverage 8 \
    --filter-total-only-biallelic-snps \
    --sync-path $SYNC \
    --out-dir $OUTPUT \
    --log-file $OUTPUT/grenedalf_fst_log.txt \
    --allow-file-overwriting \
    --threads 4