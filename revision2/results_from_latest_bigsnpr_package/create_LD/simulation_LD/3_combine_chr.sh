#!/bin/bash
#SBATCH --job-name combine
#SBATCH --mem=80G
#SBATCH --array=1-5
#SBATCH --exclude=compute-113,compute-116
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

ind1=$SLURM_ARRAY_TASK_ID

mkdir ./3_combine_chr

module load R

runr(){
    R CMD BATCH --no-save --no-restore "$1" 3_combine_chr.R ./3_combine_chr/${ind1}.Rout
}
runr "--args temp1='${ind1}'"

