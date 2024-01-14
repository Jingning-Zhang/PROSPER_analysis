#!/bin/bash
#SBATCH --job-name ldp_w
#SBATCH --mem=10G
#SBATCH --array=529-540,589-600,649-660
#SBATCH --exclude=compute-070,compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu


ind=$SLURM_ARRAY_TASK_ID

mkdir ./4_weighted

module load R

runr(){
    R CMD BATCH --no-save --no-restore "$1" 4_weighted.R ./4_weighted/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.out