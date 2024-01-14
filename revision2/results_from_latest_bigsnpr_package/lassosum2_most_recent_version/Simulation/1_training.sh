#!/usr/bin/env bash
#$ -N ls_training
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -m e
#$ -t 601-900
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
