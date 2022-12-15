#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


skipnode='!(compute-053|compute-054|compute-113)'

for setting in "5" "4" "3" "2" "1"
do

for ethnic in "AFR" "AMR" "EAS" "EUR" "SAS"
do
cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_${setting}/4_tuning_and_validation/${ethnic}
rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout
qsub -l h=${skipnode} -hold_jid three_run_${ethnic}_GA_${setting} 4_tuning_and_validation.sh
done
done
