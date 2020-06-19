#!/bin/bash

#SBATCH --job-name=extract_unaligned
#SBATCH --time=20:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --output=o_e/extract.unalign.%a.out
#SBATCH --error=o_e/extract.unalign.%a.err

declare INROOT=/data/darking1/GTEx/ncbi/version8.wgs
declare SRR_LIST=$INROOT/full.cramlist.lst

readarray -t myArray <$SRR_LIST

echo 'Downloading' ${myArray[$SLURM_ARRAY_TASK_ID]}

declare "cramfile"="${myArray[$SLURM_ARRAY_TASK_ID]}"
declare sampid=$(echo $cramfile | cut --delimiter "/" --fields 5)
declare sampid_only=$(echo $sampid | cut --delimiter "." --fields 1)

# this takes a while
ml samtools
samtools view -C -f 4 crams/$sampid > $INROOT/complete.crams/unaligned/${sampid_only}.unaligned.cram

# sbatch --array=1-452%10 extract.unaligned.sh