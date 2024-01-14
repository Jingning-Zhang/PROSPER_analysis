#!/usr/bin/env bash
#$ -N ldpredEUR
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -t 1-20
#$ -M jzhan218@jhu.edu

ind=$SGE_TASK_ID

mkdir ./4_EUR_ldpred2

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1" 4_EUR_ldpred2.R ./4_EUR_ldpred2/${ind}.Rout
}
runr "--args temp1='${ind}'"

rm *.e*
rm *.o*
