#!/bin/bash

#$ -cwd
#$ -l mem_free=13G,h_vmem=14G
#$ -e errors/
#$ -o outputs/

# This script submits R jobs to the cluster
module load perl

# please use conda_R

# module unload conda_R
# module load R/3.6.1

ID=$SGE_TASK_ID
echo ID

Rscript ./run.setup.v8.all.tissues.R ${ID}

# qsub -t 1-55 -tc 10 Rsub.1.sh

