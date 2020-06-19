#!/bin/bash

#$ -cwd
#$ -l mem_free=6G,h_vmem=7G
#$ -e errors/
#$ -o outputs/

# This script submits R jobs to the cluster
module load perl

# please do not use conda_R

module unload conda_R
module load R/3.6.1

ID=$SGE_TASK_ID
echo $ID

Rscript ./permute.GO.all.tissues.R ${ID}

# qsub -t 1-54 -tc 10 Rsub.4.sh

