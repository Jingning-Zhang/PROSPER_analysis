#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

setting='1'

skipnode='!(compute-053|compute-054|compute-113)'

for ethnic in "AFR" "AMR" "EAS" "EUR" "SAS"
do

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_${setting}/3_run_lassosum2/${ethnic}
rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout
qsub -l h=${skipnode} -hold_jid two_sumdata_${ethnic}_GA_${setting} 3_run_lassosum2.sh

done

