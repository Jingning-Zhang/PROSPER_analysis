#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


skipnode='!(compute-048|compute-053|compute-054|compute-057|compute-066|compute-068|compute-070|compute-076|compute-113)'

for setting in "1"
do

for ethnic in "AFR" "AMR" "EAS" "EUR" "SAS"
do
cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/lassosum2_sl/simulation/Setting_${setting}/${ethnic}
rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout
qsub -l h=${skipnode} simulation_1_sl.sh
done
done
