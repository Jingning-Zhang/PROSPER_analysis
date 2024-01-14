#!/bin/bash
#SBATCH --job-name ldp_training
#SBATCH --mem-per-cpu=30G --cpus-per-task=3
#SBATCH --array=437-600
#SBATCH --exclude=compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

ind=$SLURM_ARRAY_TASK_ID

mkdir ./1_training

module load R

runr(){
    R CMD BATCH --no-save --no-restore "$1" 1_training.R ./1_training/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.out
