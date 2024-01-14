#!/bin/bash
#SBATCH --job-name ls_tt
#SBATCH --mem=10G
#SBATCH --array=301-900
#SBATCH --exclude=compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

ind=$SLURM_ARRAY_TASK_ID

mkdir ./3_tuning_testing

module load R

runr(){
    R CMD BATCH --no-save --no-restore "$1" 3_tuning_testing.R ./3_tuning_testing/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.out
