#!/bin/bash
#SBATCH --job-name ls_w
#SBATCH --mem=10G
#SBATCH --array=241-720
#SBATCH --exclude=compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu


ind=$SLURM_ARRAY_TASK_ID

mkdir ./4_weighted_sl

module load R

runr(){
    R CMD BATCH --no-save --no-restore "$1" 4_weighted_sl.R ./4_weighted_sl/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.out