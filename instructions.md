# Quick cheatsheet for bash scripts
In order to execute this pipeline some basic knowledge in the use of bash scripts is expected.

## Hints to get you started

### Move between directories
Move between directories with cd (change directory) command

You can always check in which directory you are situated with the pwd (print working directory) command
```
pwd
```
### Absolute vs Relative directory names 


# How to run this pipeline
This pipeline expects a specific file structure in order to function.

...insert image ```tree -L 2```

All data files are expected to reside in the directory "data". Within the "data" directory, files are organized according to the following conventions:
Raw fastq files should be placed in the directory "untrimmed_fastq".
Reference genomes should be placed in the directory "reference".

All scripts should reside in the directory "analysis".

# Using software via Apptainer/Singularity
A container with the name poolseq_tools has already been developed with all necessary software installed.
You can launch this Singularity container by typing:
```
singularity shell poolseq_tools
```

Another way of using the container is to directly ask for a script to be run within it.
```
singularity exec poolseq_tools bash my_script.sh
```

### check quay.io to upload the image