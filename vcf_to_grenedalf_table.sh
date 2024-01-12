#!/bin/bash

# Nikos Tourvas & Lars Littmann 
# 2023.11.14
# Take a VCF produced by VarScan and output a table containing the allele frequencies for all samples per position.

# declare variables
VCF=${1}

#Take the name of the input vcf without file type suffix.
OUTPUT=${VCF/.vcf.gz/}

#Extract the number of samples found in the vcf.
N_SAMPLES=$(bcftools query -l ${VCF} | wc -l)

# We extract the following information from the VCF: 
# Chromosome, Position, Ref Allele, Alt Allele, Ref Allele read depth (per pool), Alt Allele read depth (per pool)
# bcftools query arguments (extract data from VCF)
# -H: include header
# -f: extract FORMAT columns

# After the table is produced, it is fed to awk for transformations. These transformations are:
# 1.Fix table header: (a) remove square brackets and replace first four column names 
#                     (b) replace the suffix ":RD" with ".ref.cnt" in columns which signify ref allele counts
#                     (c) replace the suffix ":AD" with ".alt.cnt" in columns which signify alt allele counts
# 2.Replace periods (.) which signify missing data with zeros (0) in the allele count columns

bcftools query -Hf '%CHROM\t%POS\t%REF\t%ALT[\t%RD][\t%AD]\n' ${VCF}\
   | awk 'BEGIN{OFS=FS="\t"} 
        NR==1 {
            $1="chrom"; $2="pos"; $3="ref"; $4="alt";
            for (i=5; i<=NF; i++){
                gsub(/\[[^]]*\]/, "", $i); 
                gsub(/:RD$/, ".ref.cnt", $i); 
                gsub(/:AD$/, ".alt.cnt", $i)
            };
        print; next} 
        NR>1 {for (i=2; i<=NF; i++) 
            gsub(/\./, NA, $i); print
        }' \
| sed -e 's/\t\t/\tNA\t/g' > ${OUTPUT}_intermediate.txt
 #Make sure that empty fields become 'NA'.

#Calculate the allele frequencies using the pairs of RD and AD fields
#The intermediate table is structured such that for each genomic position, there are a few columns with position information.
#These are followed by N_SAMPLE number of columns containing the RD values for all samples.
#Following the RD columns, there are N_SAMPLE number of columns containing the AD values for all samples in the same order as the RD values.
#Frequency values are calculated from RD and AD fields that are N_SAMPLE number of columns apart.
#This is done using a while loop

n=0 #Call a variable to store how many times the loop has run

while [ ${n} -le ${N_SAMPLES} ]; do #Run the loop as many times as there are samples in the VCF.
RD_POS=$((5+${n}));
AD_POS=$((${RD_POS}+${N_SAMPLES}));
FREQ_POS=$((${AD_POS}+${N_SAMPLES})); #Tell awk in which column to place the allele frequency it calculates.
awk -F '\t' -v awkRDpos="${RD_POS}" -v awkADpos="${AD_POS}" -v awkFREQpos="${FREQ_POS}" ' OFS = "\t" {awkFREQvalue=$awkADpos/($awkRDpos+$awkADpos); print $0,  $awkFREQpos=awkFREQvalue }' ${OUTPUT}_intermediate.txt > ${OUTPUT}_buffer.txt; #Write to an intermediate file so that awk does not try to write to the same file it is reading from
cat ${OUTPUT}_buffer.txt > ${OUTPUT}_intermediate.txt; #Replace the table with the intermediate table we just created.
echo ${n}; #print progress
n=$(( ${n}+1 )); #update how many times the loop has run
done

#Several final refinements to the table are necessary

awk -F "\t" ' OFS = "\t" {chrompos=$1"_"$2; print chrompos, $0} ' ${OUTPUT}_intermediate.txt > ${OUTPUT}_big.txt #create a column that describes the genomic position in just one cell
cut -f1,$((6+(2*${N_SAMPLES})))-$((5+(3*${N_SAMPLES}))) ${OUTPUT}_big.txt > ${OUTPUT}_small.txt #Extract only the columns that summarise the genomic position and all the allele frequencies.

head -n 1 ${OUTPUT}_big.txt | cut -f1,6-$((5+${N_SAMPLES})) | sed -e 's/P01-...-ACORN-BOKU-...-//g' | sed -e 's/.ref.cnt//g' > ${OUTPUT}_AlleleFrequencyTable.txt

sed -e 's/-nan/NA/g' ${OUTPUT}_small.txt | tail +2 >> ${OUTPUT}_AlleleFrequencyTable.txt #Make sure that NAs are noted correctly.

#clean up intermediate files
rm ${OUTPUT}_buffer.txt
rm ${OUTPUT}_intermediate.txt
