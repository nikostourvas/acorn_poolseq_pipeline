#!/bin/bash

# create output directory
mkdir -p ../results/trimmed_fastq_Batch2

# variables
IND=$1 #The individual (or pool) that will be considered in this run. GNU parallel will be used to 'assign' a pool/individual to run on.
UNTRIMMED=/mnt/AcornSeqdata/AdapterClipped_Batch2 #A directory that contains the untrimmed files that you'd like to process.
TRIMMED=/mnt/results/trimmed_fastq_Batch2 #The output directory for trimmed reads.
#ADAPTER=/mnt/data/adapters/TruSeq3-PE-2.fa #Location of a file that contains information about the adapters that were used to sequence the reads. DERELICT!!!

#The structure of the trimmomatic script is as follows (line by line). For now, the ILLUMINACLIP option is left out. 

#Call the trimomatic tool and specify that it needs to run in paired end mode (PE)
#Specify the locations of the forward and reverse untrimmed reads (seperated by a space)
#Specify the output locations and file names for the trimmed forward reads. First paired (1), and then unpaired (1un).
#Specify the output locations and file names for the trimmed reverse reads. First paired (2), and then unparied (2un).
#Invoke the option to trim Illumina adapters, followed by where to find the adapter sequences, how many mismatches are allowed to identify an adapter, in what window trimomatic should look for palindromic matches, and the threshold score for full alignments against adapters.
#Invoke the option to check for basepair read quality (Phred score) in windows. In this case, we observe windows of 4 basepairs and will only accept it if the average Phred score is above 20.
#Invoke the option to minimum length to discard reads that are deemed too short. In our case, the minimum we set is 25 reads. Any shorter increases the chances of excessive multi-mapping.

java -jar /usr/share/java/trimmomatic-0.39.jar PE \
        ${UNTRIMMED}/${IND}_R1_clipped.fastq.gz ${UNTRIMMED}/${IND}_R2_clipped.fastq.gz \
        ${TRIMMED}/${IND}_1.trim.fastq.gz ${TRIMMED}/${IND}_1un.trim.fastq.gz \
        ${TRIMMED}/${IND}_2.trim.fastq.gz ${TRIMMED}/${IND}_2un.trim.fastq.gz \
        SLIDINGWINDOW:4:20 \
        MINLEN:25

