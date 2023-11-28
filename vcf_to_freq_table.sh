#!/bin/bash

# Nikos Tourvas
# 2023.11.14
# Create a table from a varscan-produced VCF file for analysis in grenedalf.
# We extract the following information from the VCF: 
# Chromosome, Position, Ref Allele, Alt Allele, Ref Allele read depth (per pool), Alt Allele read depth (per pool)

# After the table is produced, it is fed to awk for transformations. These transformations are:
# 1.Fix table header: (a) remove square brackets and replace first four column names 
#                     (b) replace the suffix ":FREQ" with ".var.freq" in columns
# 2.Transform variant frequency percentages into decimal numbers

# declare variables
VCF=${1}
OUTPUT=${2}

# bcftools query arguments (extract data from VCF)
# -H: include header
# -f: extract FORMAT columns
bcftools query -Hf '%CHROM\t%POS\t%REF\t%ALT[\t%FREQ]\n' ${VCF} \
    | awk 'BEGIN{OFS=FS="\t"} 
        NR==1 {
            $1="chrom"; $2="pos"; $3="ref"; $4="alt"; 
            for (i=5; i<=NF; i++){ 
                gsub(/\[[^]]*\]/, "", $i); 
                gsub(/:FREQ$/, ".var.freq", $i)
            };
        print; next} 
        NR>1 {
		for (i=5; i<=NF; i++) {
                	sub(/%/, "", $i);
                	$i = $i / 100;  # convert percentage to decimal
        	}     	
        	print
        }' \
    > ${OUTPUT}.txt
