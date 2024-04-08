#!/bin/bash

#Nikos Tourvas
#31/01/2024
#This script uses bwa mem2 to align trimmed reads to a reference genome. In ACORN WP3, these reads were already deduplicated, but that is not strictly necessary for this script to run correctly.
#The input cannot be an interleaved fastq file. It requires seperate files for forward and reverse reads.
#The output of this script is a raw BAM file. In all likeliness, it will need further procesing downstream before variant calling. It is advisable to run Qualimap reports on output.
#For more information on how to run the script, see bottom of document.

# create output directory
mkdir -p /data/genetics_tmp/results/mapped_reads

# declare variables
IND=${1} #The name of the sample (or individual/IND if you will) that you are processing with this script.
REF=/data/genetics_tmp/REFERENCE/Qrob_PM1N_Organelles.fa #The path to a reference sequence. In this case the Plomion et al (2018) Quercus robur reference genome with mitochondrial and chloroplast genomes added. 
FORWARD=/data/genetics_tmp/results/fastp_dedup_trim/${IND}_1.trim.dedup.fastq.gz #The location of the fastq file containing all trimmed forward reads. Which individual is automatically specified.
REVERSE=/data/genetics_tmp/results/fastp_dedup_trim/${IND}_2.trim.dedup.fastq.gz #The location of the fastq file containing all trimmed reverse reads. Which individual is automatically specified.
RG="@RG\tID:${IND}\tPL:Illumina\tSM:${IND}" #A string that is used to specify the readgroup information in the output BAM file.
OUTPUT=/data/genetics_tmp/results/mapped_reads/${IND} #The desired location and name of the output file. Will be 'completed' in upcoming lines of the script. 
BWAMEM2=/usr/local/bin/bwa-mem2/bwa-mem2 #path to executable needs to be explicit

# Align to reference genome and export raw bam file

# BWA mem2 arguments
# -M mark shorter split hits as secondary
# -t how many cores to use PER SAMPLE for mapping

# Samtools view arguments
# -h Include header in output
# -b output a bam file

# Be mindful to account for as many cores have been assigned to bwa and samtools 
# (e.g. 2c BWA + 1c samtools = 3 cores in total)
# samtools flagstat is run afterwards 

${BWAMEM2} mem -R ${RG} -M -t 1 ${REF} ${FORWARD} ${REVERSE} 2> ${OUTPUT}.bwa-mem.err \
    | samtools view --threads 1 -h -b -o ${OUTPUT}.raw.bam 2> ${OUTPUT}.sam-view.err

# gather statistics
samtools flagstat -@ 1 ${OUTPUT}.raw.bam > ${OUTPUT}.raw.flagstat

##############################################################
#A handy way to run this script is as follows:
#Create a .txt file that contains all the sample names of samples you wish to process. These sample names have to match with your file naming. 
#Run the following parallel command. It initiates this bash script for each sample. As soon as one script finishes running, it starts the next until it has cycled through all file names in the .txt file.
#These commands are not part of the script, but can be copied into the command line to run the script. 

#parallel --verbose -j 60 \
#	'bash map.sh {}' :::: ../AcornSeqdata/samplenames.txt
