#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-066|compute-068|compute-076|compute-113)'

para='10'

M='5'

ethnic1="AFR_AMR_EAS_EUR_SAS"
ethnic2='c("AFR","AMR","EAS","EUR","SAS")'

for (( setting = 5; setting > 0; setting-- )); do
for (( rep = 1; rep < 11; rep++ )); do

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth
runr(){
    R CMD BATCH --no-save --no-restore "$1" 8_write_R_multi-lassosum_clean_score.R 8_write_R_multi-lassosum_clean_score_${M}ethnicity_${ethnic1}.Rout
}
runr "--args ethnic=${ethnic2} setting='${setting}' rep='${rep}' L='${para}' Lc='${para}' para='${para}'"

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_${setting}/${ethnic1}/8_multi-lassosum_clean_score_${para}

rm *.e*
rm *.o*
rm *.pe*
rm *.po*
#rm *.Rout
rm *_fixEUR.Rout
#qsub -l h=${skipnode} -hold_jid six_meprs-${M}ethnicity-${ethnic1}-GA_${setting}_rep_${rep}_sub rep_${rep}.sh
qsub -l h=${skipnode} -hold_jid six_meprs-${M}ethnicity-${ethnic1}-GA_${setting}_rep_${rep}_sub_fixEUR rep_${rep}_fixEUR.sh

done
done


