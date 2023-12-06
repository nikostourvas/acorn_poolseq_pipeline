#!/bin/bash

# Nikos Tourvas
# 2023.11.14
# Create a table from a varscan-produced VCF file for analysis in grenedalf.
# We extract the following information from the VCF: 
# Chromosome, Position, Ref Allele, Alt Allele, Ref Allele read depth (per pool), Alt Allele read depth (per pool)

# After the table is produced, it is fed to awk for transformations. These transformations are:
# 1.Fix table header: (a) remove square brackets and replace first four column names 
#                     (b) replace the suffix ":RD" with ".ref.cnt" in columns which signify ref allele counts
#                     (c) replace the suffix ":AD" with ".alt.cnt" in columns which signify alt allele counts
# 2.Replace periods (.) which signify missing data with zeros (0) in the allele count columns

# declare variables
VCF=${1}

OUTPUT=${VCF/.vcf/}

N_SAMPLES=$(bcftools query -l ${VCF} | wc -l)

# bcftools query arguments (extract data from VCF)
# -H: include header
# -f: extract FORMAT columns
bcftools query -Hf '%CHROM\t%POS\t%REF\t%ALT[\t%RD][\t%AD]\n' ${VCF} \
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

n=0

while [ ${n} -le ${N_SAMPLES} ]; do
RD_POS=$((5+${n}));
AD_POS=$((${RD_POS}+${N_SAMPLES}));
FREQ_POS=$((${AD_POS}+${N_SAMPLES}));
awk -F '\t' -v awkRDpos="${RD_POS}" -v awkADpos="${AD_POS}" -v awkFREQpos="${FREQ_POS}" \
' OFS = "\t" {awkFREQvalue=$awkADpos/($awkRDpos+$awkADpos); print $0,  $awkFREQpos=awkFREQvalue }' ${OUTPUT}_intermediate.txt > \
${OUTPUT}_buffer.txt;
cat ${OUTPUT}_buffer.txt > ${OUTPUT}_intermediate.txt;
echo ${n};
n=$(( ${n}+1 ));
done

awk -F "\t" ' OFS = "\t" {chrompos=$1"_"$2; print chrompos, $0} ' ${OUTPUT}_intermediate.txt > ${OUTPUT}_big.txt
cut -f1,$((6+(2*${N_SAMPLES})))-$((6+(3*${N_SAMPLES}))) ${OUTPUT}_big.txt > ${OUTPUT}_small.txt

head -n 1 ${OUTPUT}_big.txt | cut -f1,6-$((6+${N_SAMPLES})) | sed -e 's/P01-...-ACORN-BOKU-...-//g' |\
sed -e 's/.ref.cnt//g' > ${OUTPUT}_AlleleFrequencyTable.txt

sed -e 's/-nan/NA/g' ${OUTPUT}_small.txt | tail +2 >> ${OUTPUT}_AlleleFrequencyTable.txt

rm ${OUTPUT}_buffer.txt
rm ${OUTPUT}_intermediate.txt
