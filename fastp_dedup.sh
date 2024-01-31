#!/bin/bash

#31.01.24
#Nikos Tourvas & Lars Littmann
#Process adapter-clipped fastq files. The script both removes duplicate reads and trims poor-quality reads or ends of reads. 

# variables
SAMPLE=$1 #The names of directories that contain the fastq-file pairs, as provided by LGC.
INPUT_DIR=../Acorn_SeqData/AdapterClipped_Batch2 #The directory that contains all untrimmed fastq files
OUT_DIR=../results/fastp_dedup_trim #Directory in which to place all deduped and trimmed outputs.

mkdir -p ${OUT_DIR} #Make sure that the output directory exhists before running.

#from top to bottom:
#First implementation of fastp. Read the two fastq files for a particular pool/individual
#Ensure that in the first implementation there is nothing that interferes with duplicate removal. No paired read trimming
#Implement deduplication, set the accuracy to level 3 (see fastp manual), which should be good for 0.01% chance of accidental hash collision.
#Implement overlap analysis; correct any faulty basepairs in overlap regions of paired-end reads.
#Output a HTML report for deduplication. Set the title.
#Run this implementation of fastp on X threads
#Pass on the output to standardout for the next implementation of fastp (data is sent as interleaved)
#Initiate the second implementation of fastp by reading from stdin (data is received as interleaved)
#This line simulates the sliding window approach of trimomatic. Evaluate X basepairs at a time, shifting by one each time. Evaluate whether avg phred score is above y
#Remove any reads that are shorter than 50bp after trimming
#Output a HTML report for trimming. Set the title.
#Output two deduplicated and trimmed fastq files. 

fastp --in1 ${INPUT_DIR}/${SAMPLE}/${SAMPLE/Sample_/}_R1_clipped.fastq.gz --in2 ${INPUT_DIR}/${SAMPLE}/${SAMPLE/Sample_/}_R2_clipped.fastq.gz \
      --disable_adapter_trimming \
      --dedup --dup_calc_accuracy 3 \
      --correction \
      --html "${OUT_DIR}/${SAMPLE/Sample_/}_fastp_dedup.html" \
      --json "${OUT_DIR}/${SAMPLE/Sample_/}_fastp_dedup.json" \
      --thread 2 \
      --stdout |
fastp --stdin --interleaved_in \
      --disable_adapter_trimming \
      --dont_eval_duplication \
      --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
      --length_required 50 \
      --thread 2 \
      --html "${OUT_DIR}/${SAMPLE/Sample_/}_fastp_trim.html" \
      --html "${OUT_DIR}/${SAMPLE/Sample_/}_fastp_trim.json" \
      --out1 ${OUT_DIR}/${SAMPLE/Sample_/}_1.trim.dedup.fastq.gz --out2 ${OUT_DIR}/${SAMPLE/Sample_/}_2.trim.dedup.fastq.gz



      
