#!/usr/bin/env bash
#$ -N prscsx_5eth
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=100G,h='compute-119'
#$ -m e
#$ -t 1-10
#$ -M jzhan218@jhu.edu

mkdir ./prscsx_5eth

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1"  prscsx_5eth.R ./prscsx_5eth/$SGE_TASK_ID.Rout
}
runr "--args $SGE_TASK_ID"

#-l mem_free=2G,h_vmem=2G,h_fsize=100G,h='(compute-119|compute-124|compute-129)'
