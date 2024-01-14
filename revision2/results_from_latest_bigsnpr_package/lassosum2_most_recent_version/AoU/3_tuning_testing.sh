#!/usr/bin/env bash
#$ -N AoU_lassosum
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -t 1-6
#$ -M jzhan218@jhu.edu

ind=$SGE_TASK_ID

mkdir ./3_tuning_testing

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1" 3_tuning_testing.R ./3_tuning_testing/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.e*
rm *.o*
