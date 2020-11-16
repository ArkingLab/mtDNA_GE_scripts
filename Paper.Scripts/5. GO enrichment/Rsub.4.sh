#!/bin/bash

#$ -cwd
#$ -l mem_free=9G,h_vmem=10G
#$ -e errors/
#$ -o outputs/

# This script submits R jobs to the cluster
module load perl

# please do not use conda_R

module unload conda_R
module load R/3.6.1

ID=$SGE_TASK_ID
echo $ID

for i in {1..9}
do
   Rscript ./permute.GO.all.tissues_1000.R ${ID} $i
done

# qsub -t 1-54 -tc 10 Rsub.4.sh

