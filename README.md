# acorn_poolseq_pipeline

This pipeline is designed for the ingestion of short-read sequencing data and performs all the necessary steps (quality assessment, quality trimming, mapping to reference genome, variant calling) for the creation of standard VCF files which hold the information of all detected variants in the provided data set. VCF files can then be used for the exploration of adaptive genetic diversity through the outlier detection method provided by pcadapt R package.

## Software install
The required software can either be directly installed on a server, or be installed through a miniconda virtual environment. A guide to recreate a miniconda environment is provided (see miniconda_software_install.md).

## Running the pipeline
The pipeline execution is achieved from the "draft_pipeline.sh" script. Before running make sure that parameters of all nested scripts are correctly set. 

Time Benchmark: The whole pipeline runs for about 5h and 30min utilizing 30 CPU cores with peak RAM usage of around 28G, when used for the analysis of a pair of short-read pooled samples.

## Acknowledgment
This pipeline was developed with support from COST Action Conserve Plants (CA18201) during a STSM titled "Use of genomic tools for the management of forest reproductive material under climate change" under the supervision of Charalambos Neophytou at University of Natural Resources and Life Sciences (BOKU), Vienna.