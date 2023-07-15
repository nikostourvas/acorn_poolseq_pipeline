#!/bin/bash

OUTPUT=../results/frequencies

grenedalf simulate \
    --format sync \
    --random-seed 46 \
    --coverages 50,50,50,50 \
    --mutation-count 1000 \
    --length 50000 \
    --out-dir $OUTPUT \
    --compress \