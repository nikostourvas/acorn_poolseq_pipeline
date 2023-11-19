#!/bin/bash
VCF=${1} # Name of the VCF file
SUBSAMPLE_RATIO=${2:-0.01} # Default value is 0.01, can be overridden by providing a second argument

# Randomly sample the VCF to get a small working file
# This requires vcflib which is proving difficult to install correctly in our container
# bcftools view .vcf.gz | vcfrandomsample -r 0.012 > subset.vcf

# Extract Depth and Genotype Quality per pool
bcftools query -Hf '%CHROM\t%POS[\t%DP][\t%GQ]\n' ${VCF} \
    | awk 'BEGIN{OFS=FS="\t"} 
        NR==1 {
            $1="chrom"; $2="pos"; 
            for (i=3; i<=NF; i++){ 
                gsub(/\[[^]]*\]/, "", $i);
                gsub(/:DP$/, "_DP", $i);
        	    gsub(/:GQ$/, "_GQ", $i)
            };
        print; next} 
        NR>1 {
            for (i=2; i<=NF; i++) {
                gsub(/\./, "NA", $i);
            }
            print
        }' \
        | perl -ne 'print if (rand() < '${SUBSAMPLE_RATIO}' or /^chrom/)' \
    > ${VCF}_filtering_stats.tsv       