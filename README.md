# acorn_poolseq_pipeline

This pipeline is designed for the ingestion of short-read sequencing data and performs all the necessary steps (quality assessment, quality trimming, mapping to reference genome, variant calling) for the creation of standard VCF files which hold the information of all detected variants in the provided data set. VCF files can then be used for the exploration of adaptive genetic diversity through the outlier detection method provided by pcadapt R package.

## Running the pipeline
The pipeline execution is achieved from the "draft_pipeline.sh" script. This is the master script through which the remaining scripts are run. Before running make sure that parameters of all nested scripts are correctly set. 

The recommended way to run the pipeline is by using the poolseq_tools software container https://github.com/nikostourvas/poolseq_tools/tree/main