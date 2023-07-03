#!/bin/bash

#SBATCH --time=06:00:00               # Run time (days-hh:mm:ss) - (max 7days)
#SBATCH --partition=rome             # Submit to queue/partition named batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_50,TIME_LIMIT_90
#SBATCH --mail-user=tourvasn@for.auth.gr

# Load modules
module load gcc/9.2.0 singularity/3.6.4

# Load Apptainer/Singularity container
singularity exec --bind /mnt/forgenet_a/projects/ngs_training:/mnt \
        /mnt/forgenet_a/container_images/poolseq_tools_0.1.15.sif \
        /bin/sh -c 'cd /mnt/acorn_poolseq_pipeline/ && /usr/bin/time -v bash draft_pipeline.sh'