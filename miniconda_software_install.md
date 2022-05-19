# Miniconda install
first install miniconda itself
Miniconda does NOT need administrative rights to be installed.
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh
```
 
After install, we need to enable the repositories of bioinfo software. Executing these commands in the following order is important!!!
```
conda config --add channels defaults
config --add channels bioconda
conda config --add channels conda-forge
```

create new conda environment and install software
ATTN! Installing all software packages at once is the best way to resolve dependency issues.
```
conda create --name bio ea-utils fastqc multiqc trimmomatic seqtk flash bwa bowtie bedtools samtools bcftools freebayes gatk picard vcftools plink varscan primer3
```

You can start this new environment by :
```
conda activate bio
```
