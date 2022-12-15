#!/usr/bin/env bash
#$ -N LD
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

for ethnic in "AFR" "AMR" "EAS" "EUR" "SAS"
do
cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/1_prepare_LD
qsub -l h=${skipnode} ${ethnic}.sh
done
