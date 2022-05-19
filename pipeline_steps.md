# fastqc - Quality assessment
Fastqc runs in parallel by default
From the folder where the fastq or fastq.gz files are stored run
```
fastqc -t 4 -o ../../results/fastqc_untrimmed_reads/ *fastq
```

Working with the FastQC text output
```
cat */summary.txt > ~/docs/fastqc_summaries.txt

cd ~/docs/
grep FAIL fastqc_summaries.txt
```

# Trimmomatic - Quality trimming


# FLASH - How to run it in parallel for each sample
create sample list
```
for i in /home/tourvasn/ngs_training/data/raw/*_R1_001.fastq.gz; do echo $basename(${i%_R*}); done > inds
```

Finally run the **parallel_flash.sh** script with GNU parallel
```
parallel 'bash parallel_flash.sh {}' :::: inds
```

# BWA - Map to reference genome
## Reference genome
### Example download from NBCI DB
First download ncbi’s command line tool via miniconda
```
conda install -c bioconda entrez-direct 
```
Then for _Quercus robur_ 
```
esearch -db nucleotide -query “PRJEB51283” | efetch -format fasta > PRJEB51283.fasta
```

### Download latest reference genome of _Quercus robur_ from Plomion et al. 2018
```
cd ~/ngs_training/data/reference/
wget https://urgi.versailles.inra.fr/download/oak/Qrob_V2_2N.fa.gz
```

## Mapping
Go to directory where reference genome is store and type
```
bwa index Qrob_V2_N2.fa.gz
```

create sample list
```
for i in /home/tourvasn/ngs_training/data/raw/*_R1_001.fastq.gz; do echo $basename(${i%_R*}); done > inds
```

run bwa in parallel mode
```
parallel 'bash parallel_bwa.sh {}' :::: inds
```
