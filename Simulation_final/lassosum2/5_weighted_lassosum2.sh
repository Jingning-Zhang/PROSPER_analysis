#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


skipnode='!(compute-053|compute-054|compute-113)'

for setting in "5" "4" "3" "2" "1"
do
cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_${setting}/5_weighted_lassosum2
rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout
qsub -l h=${skipnode} -hold_jid four_tv_*_GA_${setting},four_tvEUR_${ethnic}_GA_${setting} 5_weighted_lassosum2.sh
done



