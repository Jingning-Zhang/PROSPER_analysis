#!/usr/bin/env bash
#$ -N AoU_lassosum
#$ -cwd
#$ -l mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -m e
#$ -t 1-6
#$ -M jzhan218@jhu.edu

ind=$SGE_TASK_ID

mkdir ./1_training

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1" 1_training.R ./1_training/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.e*
rm *.o*
