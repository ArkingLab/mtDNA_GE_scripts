#!/bin/bash

#$ -cwd
#$ -l mem_free=13G,h_vmem=14G
#$ -o outputs

# This script submits R jobs to the cluster
# please use conda_R

# module unload conda_R
# module load R/3.6.1

ID=$SGE_TASK_ID
echo ID

Rscript ./two.stage.all.blood1000.R ${ID}

echo "two stage all is finished"
# this needs to be qsubbed 100 times....how to do??? 
# you just created a for loop within the two.stage.parallel script...not the best way to solve, but it'll do.

Rscript ./two.stage.parallel.all.blood.R ${ID}

# qsub -t 55 Rsub.blood.1000.sh

