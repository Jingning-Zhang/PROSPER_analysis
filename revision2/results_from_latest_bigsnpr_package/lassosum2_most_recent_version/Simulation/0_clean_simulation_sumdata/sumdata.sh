#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G,h_fsize=100G
#$ -m e
#$ -t 1-300
#$ -M jzhan218@jhu.edu


ind1=$SGE_TASK_ID

mkdir ./sumdata

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1" sumdata.R ./sumdata/${ind1}.Rout
}
runr "--args temp1='${ind1}'"

rm *.e*
rm *.o*
