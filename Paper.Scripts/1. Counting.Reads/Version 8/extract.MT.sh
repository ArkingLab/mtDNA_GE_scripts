#!/bin/bash

#SBATCH --job-name=extract_MT
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --output=o_e/mt.extr.%a.out
#SBATCH --error=o_e/mt.extr.%a.err

declare INROOT=/data/darking1/GTEx/ncbi/version8.wgs
declare SRR_LIST=$INROOT/full.cramlist.lst

readarray -t myArray <$SRR_LIST

echo 'Extracting chrM for ' ${myArray[$SLURM_ARRAY_TASK_ID]}

declare "cramfile"="${myArray[$SLURM_ARRAY_TASK_ID]}"
declare sampid=$(echo $cramfile | cut --delimiter "/" --fields 5)
declare sampid_only=$(echo $sampid | cut --delimiter "." --fields 1)

ml samtools
samtools view -C crams/$sampid "chrM" > $INROOT/complete.crams/mito/${sampid_only}.MT.cram

# sbatch --array=1-452%20 extract.MT.sh