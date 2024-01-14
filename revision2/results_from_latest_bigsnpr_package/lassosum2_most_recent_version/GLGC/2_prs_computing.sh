#!/bin/bash
#SBATCH --job-name ls_glgc_computing
#SBATCH --mem=3G
#SBATCH --array=1-20
#SBATCH --exclude=compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

ind=$SLURM_ARRAY_TASK_ID

mkdir ./2_prs_computing

module load R

runr(){
    R CMD BATCH --no-save --no-restore "$1" 2_prs_computing.R ./2_prs_computing/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.out
