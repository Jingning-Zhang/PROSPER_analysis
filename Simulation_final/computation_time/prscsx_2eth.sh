#!/usr/bin/env bash
#$ -N prscsx_2eth
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=100G,h='compute-119'
#$ -m e
#$ -t 1-10
#$ -M jzhan218@jhu.edu

mkdir ./prscsx_2eth

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore "$1"  prscsx_2eth.R ./prscsx_2eth/$SGE_TASK_ID.Rout
}
runr "--args $SGE_TASK_ID"

