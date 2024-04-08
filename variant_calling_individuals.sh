#Lars Littmann
#2023.08.23
#Call variants (SNPs and INDELs) from BAMs that contain INDIVIDUAL sequencing data (not pool sequenced).
#For the ACORN project, these individuals are 2x20 Quercus robur trees belonging to populations 305 and 312
#This script is optimised for working on arbitrary, smaller subsections of a genome (Chunks). See the Simplified_Chunks.sh script.
#https://github.com/nikostourvas/acorn_poolseq_pipeline/blob/singularity/Simplified_Chunks.sh
#This variant calling script does not output one complete vcf, but one vcf per specified chunk. 
#Though this script is very useful for massive parallelisation, it requires a bit more manual manipulation between runs.
#Some variables need to be changed within the script before initiating a new variant calling event. These variables have been marked UPDATE BEFORE RUNNING
#A parallel command that should be run in the command line is included at the bottom of this document. 

BAM_LIST_IND=/data/genetics_tmp/VCF_Qpubescens_Pools/Qpubescens_poolBams_sorted_FastStorage.txt #For downstream convenience, we recommend that this list of bams is sorted in a logical order.
OUTDIR=/data/genetics_tmp/VCF_Qpubescens_Pools/ #Specify an output directory. UPDATE BEFORE RUNNING
REF=/mnt/reference/Qrob_PM1N.fa #Location of the reference genome that raw reads were aligned to.
CHUNK=$1 #Chunk that the current instance of the script is running on. Specified in command line. 
CHUNK_SHORT=$(basename ${CHUNK/.bed/})
THREADS=1

#Make sure the output directory is created.
mkdir -p ${OUTDIR}

#Create a .txt file that contains the file names (without directory information or suffixes) found in the BAM_LIST.
#First, create a directory to store the sample naming lists in.

mkdir -p ${OUTDIR}/SampleNamingFiles

#Then, remove any possible older versions of this file
rm ${OUTDIR}/SampleNamingFiles/SampleNaming_VCF_IND_${CHUNK_SHORT}.txt

while read line; do 
SAMPLE_NAME=$(basename ${line} | cut -d "." -f 1);
printf "%s\n" "${SAMPLE_NAME}" >> /mnt/results/SampleNaming_VCF_IND_${CHUNK_SHORT}.txt;
done < ${BAM_LIST_IND}

# Create a mpileup file for each genomic region and call snps & indels together
# Input: (i) Filtered BAM files, (ii) indexed reference genome
# Output: (i) compressed genomic region VCF files for SNPs,
#         (ii) compressed genomic region VCF files for INDELs

# Multiple instances of this script can be run concurrently with GNU parallel.
# The mpileup | varscan | bgzip pipeline requires 2 CPU threads to run
# efficiently. For this reason we assign 2 CPU threads to the downstream
# bcftools steps.

# samtools mpileup arguments:
# -B: disable Base Alignment Quality (BAQ) adjustment as recommended by VarScan
# publication (Kobolt et al., 2013)
# -q: minimum mapping quality
# -r: region where variant calling is performed
# -l: Specific chunk within region where variant calling is performed
# -f: reference genome file

# varscan arguments:
# --vcf-sample-list: a list of sample names in the same order as BAM files,
# one per line
# --min-coverage: Minimum read depth at a position to make a call
# --min-var-freq: Minimum variant allele frequency threshold
# --min-reads2: Minimum supporting reads at a position to call variants
# --min-freq-for-hom: Minimum frequency to call homozygote
# --p-value: P-value threshold for variant calling. It is recommended to use
# less stringent values for INDELs than SNPs. Here we call INDELs & SNPs
# together so we use p=0.1 for both. We will filter SNPs at a later
# stage though with the script "snp_indel_rm.sh".
# --output-vcf: Set to 1, to produce VCF file instead of table of alleles

samtools mpileup -B -q 20 -l ${CHUNK} -f ${REF} -b ${BAM_LIST_IND} -o ${OUTDIR}/${CHUNK_SHORT}_samtools.mpileup \
    2> ${OUTDIR}/${CHUNK_SHORT}_ind.mpileup.err &&
    
java -jar /usr/share/java/varscan.jar mpileup2cns ${OUTDIR}/${CHUNK_SHORT}_samtools.mpileup\
    --min-coverage 5 \
    --min-var-freq 0.025 \
    --min-reads2 1 \
    --min-freq-for-hom 0.85 \
    --vcf-sample-list /mnt/results/SampleNaming_VCF_IND_${CHUNK_SHORT}.txt \
    --p-value 0.1 \
    --output-vcf 1 \
    2> ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.err \
    | bgzip --compress-level -1 2> ${OUTDIR}/${CHUNK_SHORT}_Gzipping.err\
        > ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.vcf.gz &&

# index vcfs
bcftools index ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.vcf.gz \
    --threads ${THREADS} \
    2> ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.snpindel.vcf.err

# extract snps and save them in a separate compressed VCF
bcftools view -v snps --threads ${THREADS} \
    ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.vcf.gz \
    -Oz -o ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snp.vcf.gz \
    2> ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_view.snp.vcf.err \

# extract indels and save them in a separate compressed VCF
bcftools view -v indels --threads ${THREADS} \
    ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.vcf.gz \
    -Oz -o ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.indel.vcf.gz \
    2> ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools.indel.vcf.err

# remove redundant files for storage efficiency
rm ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.vcf.gz \
   ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.vcf.gz.csi

# some scaffolds will have no snps/indels
# delete the empty vcf files originating from these scaffolds
#find ${OUTDIR}/${CHUNK_SHORT}*.gz -maxdepth 1 -type f -empty -print -delete
#Defunct, but not tested since.

# index newly created VCFs
bcftools index ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snp.vcf.gz \
    --threads ${THREADS} \
    2> ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.snp.vcf.err

bcftools index ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.indel.vcf.gz \
    --threads ${THREADS} \
    2> ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.indel.vcf.err &&

#The following commands are used to output a single, well-formatted log file for each of the chunks, to avoid overcrowding directories.
echo -e "All log files for ${CHUNK_SHORT/.bed/}\n\n#####\n\nSamtools mpileup\n\n">${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.mpileup.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nvarScan.snpindel\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.varScan.snpindel.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools index\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.snpindel.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools view snp\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_view.snp.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools view indels\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools.indel.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools index snp\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.snp.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
echo -e "\n\n#####\n\nbcftools index snp\n\n" >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log
cat ${OUTDIR}/${CHUNK_SHORT}_ind.bcftools_index.indel.vcf.err >> ${OUTDIR}/AllLogFiles_${CHUNK_SHORT}.log &&

#Remove intermediate files to avoid overcrowding directories.
rm ${OUTDIR}/${CHUNK_SHORT}*.err 
rm ${OUTDIR}/${CHUNK_SHORT}_samtools.mpileup

########################################################################
#A handy way to run this script is as follows:
#Create a .txt file that contains all the sample names of samples you wish to process. These sample names have to match with your file naming. 
#Run the following parallel command. It initiates this bash script for each sample. As soon as one script finishes running, it starts the next until it has cycled through all file names in the .txt file.
#These commands are not part of the script, but can be copied into the command line to run the script. 

#parallel --verbose -j 80 \
#	'bash /mnt/acorn_poolseq_pipeline/variant_calling_individuals.sh {}' :::: /data/genetics_tmp/REFERENCE/ChunkFiles/Locations_Of_Chunk_Beds.txt
