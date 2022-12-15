#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

para='10'

M='5'

ethnic1="AFR_AMR_EAS_EUR_SAS"
ethnic2='c("AFR","AMR","EAS","EUR","SAS")'

for (( setting = 5; setting > 0; setting-- )); do
  for (( rep = 1; rep < 11; rep++ )); do


cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth
runr(){
    R CMD BATCH --no-save --no-restore "$1" 2_write_R_multi-lassosum.R 2_write_R_multi-lassosum_setting${setting}_${M}ethnicity_${ethnic1}_rep_${rep}.Rout
}
runr "--args ethnic=${ethnic2} setting='${setting}' rep='${rep}' L='${para}' Lc='${para}' para='${para}'"

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_${setting}/${ethnic1}/2_multi-lassosum_by_chr_${para}

rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout

qsub -l h=${skipnode} -hold_jid five_mdata-${M}ethnicity-${ethnic1}-GA_${setting}_rep_${rep} rep_${rep}.sh

done
done

