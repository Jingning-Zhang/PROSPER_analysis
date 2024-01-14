#!/bin/bash
#SBATCH --job-name ldp_aou_eur
#SBATCH --mem=3G
#SBATCH --array=1-6
#SBATCH --exclude=compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

ind=$SLURM_ARRAY_TASK_ID

mkdir ./4_EUR_ldpred2

module load R

runr(){
    R CMD BATCH --no-save --no-restore "$1" 4_EUR_ldpred2.R ./4_EUR_ldpred2/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.out
