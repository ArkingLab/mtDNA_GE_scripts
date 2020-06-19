#!/bin/bash -l

#SBATCH --job-name=counting_files
#SBATCH --time=72:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks-per-node=8
#SBATCH --output=fail_%a.out
#SBATCH --error=fail_e_%a.err

declare INROOT=/data/darking1/GTEx/ncbi/version8.wgs
declare SRR_LIST=$INROOT/full.cramlist.lst

readarray -t myArray <$SRR_LIST

ml samtools
i=$SLURM_ARRAY_TASK_ID
declare "cramfile"="${myArray[${i}]}"
declare sampid=$(echo $cramfile | cut --delimiter "/" --fields 5)

count=$(samtools view -c /work-zfs/darking1/resources/GTEx/ncbi/version8.wgs/complete.crams/crams/${sampid_only}.cram)
echo "${count},${sampid_only}" >> read.info/more.reads2/${sampid_only}.reads.csv

# sbatch --array=1-452%20 count.full.reads.missed.sh



