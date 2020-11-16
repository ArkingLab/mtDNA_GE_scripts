#!/bin/bash

#$ -cwd
#$ -l mem_free=13G,h_vmem=14G

# This script submits R jobs to the cluster
# please use conda_R

module load conda_R/3.6.x
# module load R/3.6.1

ID=$SGE_TASK_ID
echo $ID

Rscript ./two.stage.all.R ${ID}

echo "two stage all is finished"
# this needs to be qsubbed 1000 times....how to do??? 
# you just created a for loop within the two.stage.parallel script...not the best way to solve, but it'll do.

Rscript ./two.stage.parallel.all.1000.R ${ID}

# qsub -t 1:55 -tc 10 Rsub.all.tissues.1000.sh

