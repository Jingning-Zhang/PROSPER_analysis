#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

M='2'

for (( rep = 1; rep < 11; rep++ )); do
  for (( setting = 1; setting < 6; setting++ )); do

for ethnic in "AFR" "AMR" "EAS" "SAS"
do

    ethnic1=${ethnic}_EUR

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum
runr(){
    R CMD BATCH --no-save --no-restore "$1" 1_write_R_standard_data.R 1_write_R_standard_data_setting${setting}_${M}ethnicity_${ethnic}_rep_${rep}.Rout
}
runr "--args ethnic=c('${ethnic}','EUR') setting='${setting}' rep='${rep}'"

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/Setting_${setting}/${ethnic1}/1_standard_data
rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout
qsub -l h=${skipnode} -hold_jid four_tv_*_GA_${setting},ld-* rep_${rep}.sh

done

done

done



