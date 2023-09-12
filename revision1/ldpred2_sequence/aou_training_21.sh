#!/usr/bin/env bash
#$ -N aou_21
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu


ind1='2'
ind2='1'

mkdir ./aou_training

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1" aou_training.R ./aou_training/${ind1}${ind2}-aou_training_chr$SGE_TASK_ID.Rout
}
runr "--args temp1='${ind1}' temp2='${ind2}' chr='$SGE_TASK_ID'"
