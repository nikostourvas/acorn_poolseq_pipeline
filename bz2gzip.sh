#!/bin/bash

# Nikos Tourvas
# Convert FASTQ files compressed as bz2 to gz. The original files are retained.
# The script expects that each FASTQ file is found inside a separate subdirectory.

INPUT=../data/untrimmed_fastq_ind/AdapterClipped #Path to the 'mother directory' that contains all subdirectories that contain FASTQ files.
PATTERN=s*/*.bz2 #A pattern that can be used to find all FASTQ files, regardless of them being in subdirectories.

cd ${INPUT} #Navigate to the 'mother directory' that contains all subdirectories that contain FASTQ files.
find ${PATTERN} -type f > indfile.txt #Create a textfile in the 'mother directory'. Each line in the textfile is a unique path to one of the FASTQ files we wish to process.

#The following parallel command loops through the textfile 'infile.txt' and opens the specified .bz2-zipped file defined by a line, pipes it as stdout to gzip, which in turn re-zips it as a .gz file.
parallel --dry-run --verbose -j 18 \
    'bzcat {} | gzip -c > {.}.gz' :::: indfile.txt

rm indfile.txt #Remove the file 'indfile.txt' after the script has finished running to keep the directory clutter-free. 
