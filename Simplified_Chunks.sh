#!/bin/bash

#Lars Littmann 
#2023.09.06
#Divide a reference genome up into similar sized, arbitrary 'chunks' in order to make more efficient use of paralellisation.
#The actual reference genome is not split up. Instead, this script generates a list of regions that can be used by variant calling scripts.


#Inputs: a fasta index file (.fa.fai) containing all chromosomes and scaffolds, and the desired length in bp of every chunk.

INDEX=$1 #Fasta index file (.fa.fai) of the reference genome
TARGET_SIZE=$2 #The size in bp that every chunk should roughly be.

#Create a good place to output all the .bed files

mkdir -p ../reference/ChunkFiles

rm ../reference/ChunkFiles/*unk* #Remove anything chunk-related that is already in the output directory.
#This step is necessary because this script CONCETENATES output onto a file, instead of OVERWRITING anything. 

OUTDIR=../reference/ChunkFiles

#Convert the .fai into a .bed, so we have all of our genomic coordinates.

awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${INDEX} > ${INDEX/.fa.fai/.bed} # Take info about contig name and length from .fa.fai index file and use it to construct a .bed file

bedtools makewindows -b ${INDEX/.fa.fai/.bed} -w ${TARGET_SIZE} > ${INDEX/.fa.fai/_windowed.bed} # Divide any large chromosomes or scaffolds up into chunks that match the target size.

TICKER=1;
while read line; do
printf "%s\t%s\t%s\n" ${line} > ${OUTDIR}/chunk${TICKER}.bed;
TICKER=$((${TICKER}+1))
done < ${INDEX/.fa.fai/_windowed.bed} #Use the windowed .bed as input for this loop.

#Store the paths to all of the chunk files in a convenient file that can be used later by gnu parellel when running the genotyping. 
realpath ${OUTDIR}/chunk*.bed > ${OUTDIR}/Locations_Of_Chunk_Beds.txt

#Remove some clutter files that are no longer needed.
rm ${INDEX/.fa.fai/_windowed.bed}
rm ${INDEX/.fa.fai/.bed}
