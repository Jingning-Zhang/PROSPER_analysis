#!/usr/bin/env bash
#$ -N compute_LD_4
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -m e
#$ -t 1-12
#$ -M jzhan218@jhu.edu


ind1='4'
ind2=$SGE_TASK_ID

mkdir ./2_compute_LD

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1" 2_compute_LD.R ./2_compute_LD/${ind1}-compute_LD_chr$SGE_TASK_ID.Rout
}
runr "--args temp1='${ind1}' temp2='${ind2}'"
