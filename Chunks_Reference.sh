#!/bin/bash

#Lars Littmann 
#2023.08.21
#Divide a reference genome up into similar sized, arbitrary 'chunks' in order to make more efficient use of paralellisation.
#The actual reference genome is not split up. Instead, this script generates multiple lists of genomic coordinates that can be used as input for downstream tools. 
#Each chunk is given a number. VCFs can at a later point be merged based on this number.

#Inputs: a fasta index file (.fa.fai) containing all chromosomes and scaffolds, and the desired length in bp of every chunk.

INDEX=$1 #Fasta index file (.fa.fai) of the reference genome
TARGET_SIZE=$2 #The size in bp that every chunk should roughly be.

#Create a good place to output all the .bed files

mkdir -p /home/geneticsShare/reference/ChunkFiles

OUTDIR=/home/geneticsShare/reference/ChunkFiles

#Convert the .fai into a .bed, so we have all of our genomic coordinates.

awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${INDEX} > ${OUTDIR}/${INDEX/.fa.fai/.bed} # Take info about contig name and length from .fa.fai index file and use it to construct a .bed file

bedtools makewindows -b ${OUTDIR}/${INDEX/.fa.fai/.bed} -w ${TARGET_SIZE} > ${OUTDIR}/${INDEX/.fa.fai/_windowed.bed} # Divide any large chromosomes or scaffolds up into chunks that match the target size.

#The following chunk of code is adapted from an older script of mine. I simply plug in the target size. This code spits out a bunch of numbered .bed files. 
#These bed files can be used downstream by samtools mpileup as input for the -l flag, which tells mpileup to only operate for a specific region.
TICKER1=1; #Counts up from 1 to generate numbers for the file names. The region files will be in the same order as they are in the reference genome.
TICKER2=0; #An 'empty' counter that keeps track of how many basepairs are in a region file.
MAXTICKER2=${TARGET_SIZE};
while read line; do 
if [ ${TICKER2} -lt ${MAXTICKER2} ]; #Check whether a region file has been 'filled up' to the target size.
    then printf ${line}'%s\n' >> ${OUTDIR}/chunk${TICKER1}.bed; #If there is still space for more basepairs in the region file, add the next fragment, which is in practical terms a line from the input .bed file.
    ADDTICKER=$(echo ${line} | awk '{print $3;}'); #Determine how many basepairs were added to the region file that we are currently 'filling up' in the last line.
    TICKER2=$((${TICKER2}+${ADDTICKER})); #Add how many basepairs were addded in this iteration to how many basepairs were already present.
    else TICKER2=0; TICKER1=$((${TICKER1}+1)); #If the file is found to be full enough already, reset the counter that keeps track of how many basepairs are in a region file, and update the counter that names file by adding 1.
    printf ${line}'%s\n' >> ${OUTDIR}/chunk${TICKER1}.bed; #Instead of adding the next line of the .bed file to the 'full' region file, add it to an empty one.
    ADDTICKER=$(echo ${line} | awk '{print $3;}'); #Determine how many basepairs have been added to this new file
    TICKER2=$((${TICKER2}+${ADDTICKER})); #Update the counter that keeps track of how many basepairs were added to this new file.
fi;
done <${OUTDIR}/${INDEX/.fa.fai/_windowed.bed} #Use the windowed .bed as input for this loop.

#Store the paths to all of the chunk files in a convenient file that can be used later by gnu parellel when running the genotyping. 
realpath ${OUTDIR}/chunk*.bed > ${OUTDIR}/Locations_Of_Chunk_Beds.txt

#Remove some clutter files that are no longer needed.
rm ${INDEX/.fa.fai/_windowed.bed}
rm ${INDEX/.fa.fai/.bed}
