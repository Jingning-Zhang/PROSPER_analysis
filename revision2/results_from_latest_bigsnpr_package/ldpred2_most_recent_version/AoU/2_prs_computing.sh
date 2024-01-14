#!/usr/bin/env bash
#$ -N AoU_ldpred
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -t 1-6
#$ -M jzhan218@jhu.edu

ind=$SGE_TASK_ID

mkdir ./2_prs_computing

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1" 2_prs_computing.R ./2_prs_computing/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.e*
rm *.o*
