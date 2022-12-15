#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-054|compute-053)'

setting='3'

para='5'

rep='1'
M='2'


for ethnic in "AFR" "AMR" "EAS" "SAS"
do

  ethnic1=${ethnic}_EUR

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth
runr(){
    R CMD BATCH --no-save --no-restore "$1" 2_write_R_multi-lassosum.R 2_write_R_multi-lassosum_setting${setting}_${M}ethnicity_${ethnic}_rep_${rep}.Rout
}
runr "--args ethnic=c('${ethnic}','EUR') setting='${setting}' rep='${rep}' L='${para}' Lc='${para}' para='${para}'"

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_${setting}/${ethnic1}/2_multi-lassosum_by_chr_${para}

rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout

qsub -l h=${skipnode} -hold_jid five_mdata-${M}ethnicity-${ethnic1}-GA_${setting}_rep_${rep} rep_${rep}.sh

done

