#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -t 1-5
#$ -m e
#$ -M jzhan218@jhu.edu

readarray -t a < /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/ancestry
readarray -t b < /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/trait

for i in {1..5}
do
  cd /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/${a[$(($SGE_TASK_ID-1))]}
  gunzip < ${b[$i-1]}.gz > ${b[$i-1]}
done
