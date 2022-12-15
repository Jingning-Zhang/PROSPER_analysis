#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'


M='5'

ethnic2='c("AFR","AMR","EAS","EUR","SAS")'
ethnic1="AFR_AMR_EAS_EUR_SAS"

for (( rep = 1; rep < 11; rep++ )); do
  for (( setting = 1; setting < 6; setting++ )); do

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum
runr(){
    R CMD BATCH --no-save --no-restore "$1" 1_write_R_standard_data.R 1_write_R_standard_data_setting${setting}_${M}ethnicity_${ethnic1}_rep_${rep}.Rout
}
runr "--args ethnic=${ethnic2} setting='${setting}' rep='${rep}'"

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/Setting_${setting}/${ethnic1}/1_standard_data
rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout
qsub -l h=${skipnode} -hold_jid four_tv_*_GA_${setting},ld-* rep_${rep}.sh

  done

done

