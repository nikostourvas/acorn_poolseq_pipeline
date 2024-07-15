#!/bin/bash

#Lars Littmann
#22.06.2024
#Take an AF table, Thinned AF table, and parameters for a full LFMM run that loops through multiple Ks and environmental factors. 
#FIRST Impute the AF table and the thinned AF table
#SECOND Split up the AF table to enable parallel processing
#THIRD Launch parallel lfmm jobs with all the specified parameters


###FIRST STEP###
#The imputation-inator
#Imputes missing values in an allele frequency table by calculating the mean of the allele frequency values that are there.
#Allele frequency tables are created using either Make_Frequency_Table.sh or Make_Thinned_AFtable.sh
#For the calculation of the mean, NA fields are fully ignored.
#Missingness per SNP should be kept below 10%, otherwise this imputation method is not reliable.

# declare variables
INPUT_AF_TABLE=$1 #The complete allele frequency table as produced by the script Make_Frequency_Table.sh. No need to perform any imputation etc. beforehand. Needs to be placed in ./dat
INPUT_THINNED_AF_TABLE=${INPUT_AF_TABLE/AlleleFrequencyTable.txt/Thinned_500bp_AlleleFrequencyTable.txt}  #The thinned allele frequency table as produced by Make_Thinned_AFtable.sh. 
                                                                                                          #No need to perform any imputation etc. beforehand. Needs to be placed in ./dat

ENV_DATA=$2 #A csv that contains all the environmental data (in columns) for each of the populations (in rows). Population names need to correspond with the genetic data.
SET_K=$3 #The K that you wish the script to analyse. The script will only run for one instance of K.
POPULATIONS=$4 #The 3-character code that specifies which (sub)set of ACORN populations is used.
ENV_VARIABLES=$5 #A list seperated by commas WITHOUT SPACES of the names of all the environmental factors you wish to include in the analysis.

#Create a good filename for the imputed table outputs. Will be stored as intermediate and can be used by other analyses that require imputed datasets.
IMPUTED_TABLE=$(basename ${INPUT_AF_TABLE/.txt/_Imputed_by_Mean.txt})
IMPUTED_THINNED_TABLE=$(basename ${INPUT_THINNED_AF_TABLE/.txt/_Imputed_by_Mean.txt})

#Create an output directory for this particular dataset's intermediate files
OUTPUT_DIR=$(dirname ${INPUT_AF_TABLE})/${POPULATIONS}_IntermediateFiles_LFMM
mkdir -p ${OUTPUT_DIR}

#The first line creates a new column that contains the mean of all values. 
#The second line replaces 'NA' values with the mean which is contained in the last column.
#At the end of the second line, the column containing the mean values is removed. 
#The fourth line removes any lines in which there is no variance. Required for LFMM to work.
#This takes care of an edge case where only one population has a valid allele frequency and the rest of the populations registers as NA
awk -F "\t" ' OFS="\t" {sum = 0; j = 1; MEANpos = NF+1; for (i = 2; i <= NF; i++) if ($i=='NA') {j++} else (sum+=$i); sum /= (NF-j); $MEANpos=sum; print $0 }' ${INPUT_AF_TABLE} | \
awk -F "\t" ' OFS="\t" {sum = 0; MEANpos=NF; VARpos=NF+1; for (i=2; i<=NF-1; i++) sum+=($1-$MEANpos)^2; $VARpos=sum; print $0 }' | \
awk -F "\t" ' OFS="\t" {MEANpos= NF-1; VARpos=NF; for (i=2; i <= NF; i++) if ($i=="NA") {$i=$MEANpos}; if ($VARpos!=0) {print $0}} ' | \
awk -F "\t" ' OFS="\t" {NF-=2}1' > ${OUTPUT_DIR}/${IMPUTED_TABLE}

#Repeat for the Thinned table
awk -F "\t" ' OFS="\t" {sum = 0; j = 1; MEANpos = NF+1; for (i = 2; i <= NF; i++) if ($i=='NA') {j++} else (sum+=$i); sum /= (NF-j); $MEANpos=sum; print $0 }' ${INPUT_THINNED_AF_TABLE} | \
awk -F "\t" ' OFS="\t" {sum = 0; MEANpos=NF; VARpos=NF+1; for (i=2; i<=NF-1; i++) sum+=($1-$MEANpos)^2; $VARpos=sum; print $0 }' | \
awk -F "\t" ' OFS="\t" {MEANpos= NF-1; VARpos=NF; for (i=2; i <= NF; i++) if ($i=="NA") {$i=$MEANpos}; if ($VARpos!=0) {print $0}} ' | \
awk -F "\t" ' OFS="\t" {NF-=2}1' > ${OUTPUT_DIR}/${IMPUTED_THINNED_TABLE}

###SECOND STEP###
#The splitter-upper
#Splits up the large, unthinned AFtable into 50 roughly equal-size tables. The thinned dataset is small enough as it is and should not be split.

HEADER=$(head -n 1 ${OUTPUT_DIR}/${IMPUTED_TABLE}) #Store the head line of the table in a string for later.
tail -n +2 ${OUTPUT_DIR}/${IMPUTED_TABLE} > body.txt #Create a file that contains the whole table EXCEPT for the header line.
split -n l/50 --numeric-suffixes=01 --additional-suffix .txt body.txt ${OUTPUT_DIR}/${IMPUTED_TABLE/.txt/_Chunk} #Split up the header-less table
sed -i "1i ${HEADER}" ${OUTPUT_DIR}/${IMPUTED_TABLE/.txt/_Chunk??.txt} #Insert the header line at the top of each split table
rm body.txt #Remove the header-less table that we created.

###THIRD STEP###
#The R initiation-station
#Start by creating a file that contains all the parameters for each job. Then use this file to instruct GNU parallel.

realpath ${OUTPUT_DIR}/${IMPUTED_TABLE/.txt/_Chunk??.txt} > LFMM_Full_Analysis_${POPULATIONS}_INTERMEDIATE.txt #Use the path of all 50 split genetic datasets as a scaffold to create the parameters file.

#The following awk command adds all the other parameters in the order that the next R-script expects. Separated by a space (except for the environmental variables, which are given as a comma-seperated list).
awk -F " " -v awk_working_directory="${PWD}" -v awk_thinned_dataset="${IMPUTED_THINNED_TABLE}" \
-v awk_environmental_data="${ENV_DATA}" -v awk_set_k="${SET_K}" -v awk_populations="${POPULATIONS}" -v awk_environmental_factors=${ENV_VARIABLES} \
' OFS=" " {print awk_working_directory, $0, awk_thinned_dataset, awk_environmental_data, awk_set_k, awk_populations, NR, awk_environmental_factors}' \
LFMM_Full_Analysis_${POPULATIONS}_INTERMEDIATE.txt > ./par/LFMM_Full_Analysis_Part1_Parameters_${POPULATIONS}.txt

rm LFMM_Full_Analysis_${POPULATIONS}_INTERMEDIATE.txt #Get rid of unneeded intermediate file.

###BONUS STEP###
#Create a file that can be used to launch the next Rscript as well. The user can use this file to launch the second R-script with a gnu parallel command.

echo ${ENV_VARIABLES} | tr -s ',' '\n' > LFMM_Full_Analysis_Part2_${POPULATIONS}_INTERMEDIATE.txt #This time, use the list of environmental variables as a scaffold for the file. The comma-separated variable names are properly seperated. 

#The following awk command adds all the other parameters in the order that the next R-script expects. Separated by a space.
awk -F " " -v awk_working_directory="${PWD}" -v awk_populations="${POPULATIONS}" -v awk_set_k="${SET_K}" \
'OFS=" " {print awk_working_directory, awk_populations, awk_set_k, $0}' LFMM_Full_Analysis_Part2_${POPULATIONS}_INTERMEDIATE.txt > ./par/LFMM_Full_Analysis_Part2_Parameters_${POPULATIONS}.txt

rm LFMM_Full_Analysis_Part2_${POPULATIONS}_INTERMEDIATE.txt

###THE GRANDE FINALE###
#Actually launching the first round of lfmm jobs

parallel --verbose -j 50 'Rscript /home/geneticsShare/acorn_poolseq_pipeline/LFMM_Full_analysis_Part1.R {}' :::: ./par/LFMM_Full_Analysis_Part1_Parameters_${POPULATIONS}.txt

