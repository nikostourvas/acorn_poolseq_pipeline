#!/bin/bash

# declare variables
TABLE=/data/genetics_tmp/VCF_AllPools_3Species/VCF_AllPools_3SpeciesIndelFilteredSNPs_Biallelic_grenedalf.txt
POOLSIZES=/data/genetics_tmp/variant_calling_tmp_storage_all_pools/Allpools_3species_poolsizes.txt
REFERENCE=/mnt/reference/Qrob_PM1N.fa
OUTDIR=/data/genetics_tmp/grenedalf_Allpools_3species
THREADS=12

mkdir -p ${OUTDIR} 

grenedalf fst \
    --window-type genome \
    --method unbiased-nei \
    --pool-sizes $OUTDIR/$POOLSIZES \
    --filter-sample-min-count 8 \
    --filter-sample-max-coverage 500 \
    --filter-sample-min-coverage 40 \
    --filter-total-only-biallelic-snps \
    --reference-genome-fasta-file $REFERENCE\
    --vcf-path $TABLE \
    --out-dir $OUTDIR \
    --log-file $OUTDIR/grenedalf_fst_genome_log.txt \
    --allow-file-overwriting \
    --threads $THREADS