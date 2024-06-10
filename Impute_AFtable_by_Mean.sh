#!/bin/bash

# Lars Littmann 
# 2024.06.10

#The imputation-inator
#Imputes missing values in an allele frequency table by calculating the mean of the allele frequency values that are there.
#Allele frequency tables are created using either Make_Frequency_Table.sh or Make_Thinned_AFtable.sh
#For the calculation of the mean, NA fields are fully ignored.
#Missingness per SNP should be kept below 10%, otherwise this imputation method is not reliable.

# declare variables
INPUT_AF_TABLE=$1

#Create a good filename for the output.
OUTPUT_TABLE=${INPUT_AF_TABLE/.txt/_Imputed_by_Mean.txt}

#The first line creates a new column that contains the mean of all values. 
#The second line replaces 'NA' values with the mean which is contained in the last column.
#At the end of the second line, the column containing the mean values is removed. 
awk -F "\t" ' OFS="\t" {sum = 0; j = 1; MEANpos = NF+1; for (i = 2; i <= NF; i++) if ($i=='NA') {j++} else (sum+=$i); sum /= (NF-j); $MEANpos=sum; print $0 }' ${INPUT_AF_TABLE} | \
awk -F "\t" ' OFS="\t" {MEANpos= NF; for (i=2; i <= NF; i++) if ($i=="NA") {$i=$MEANpos}; $(NF--); print $0} ' > ${OUTPUT_TABLE}
