#!/bin/bash

#$ -cwd
#$ -l mem_free=13G,h_vmem=14G
#$ -e errors/
#$ -o outputs/

# This script submits R jobs to the cluster
module load perl

# module unload conda_R
# module load R
# ID=$1
# echo ${ID}
# R CMD BATCH ./${ID}
ID=$SGE_TASK_ID
echo $ID
# Rscript ./run.setup.R ${ID}
Rscript ./run.test.all.tissues.R ${ID}

# qsub -t 1-55 -tc 10 Rsub.2.sh

