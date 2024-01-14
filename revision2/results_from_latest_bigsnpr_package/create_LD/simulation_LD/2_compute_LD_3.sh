#!/bin/bash
#SBATCH --partition=shared
#SBATCH --job-name compute_LD_3
#SBATCH --mem=20G
#SBATCH --array=1-22
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

ind1='3'
ind2=$SLURM_ARRAY_TASK_ID

mkdir ./2_compute_LD

module load R

runr(){
    R CMD BATCH --no-save --no-restore "$1" 2_compute_LD.R ./2_compute_LD/${ind1}-compute_LD_chr$SLURM_ARRAY_TASK_ID.Rout
}
runr "--args temp1='${ind1}' temp2='${ind2}'"
