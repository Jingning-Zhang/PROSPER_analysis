#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-048|compute-053|compute-054|compute-066|compute-068|compute-070|compute-076|compute-113)'

para='5'
M='5'

ethnic1='AFR_AMR_EAS_SAS_EUR'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/23andme/multi-lassosum_train_from_single_eth
runr(){
    R CMD BATCH --no-save --no-restore "$1" 6_write_R_clean_score_send_to_23andme.R 6_write_R_clean_score_send_to_23andme_${M}ethnicity_${ethnic1}.Rout
}
runr "--args ethnic=c('AFR','AMR','EAS','SAS','EUR') L='${para}' Lc='${para}' para='${para}'"


cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/23andme/multi-lassosum_train_from_single_eth/${ethnic1}/6_clean_score_send_to_23andme_${para}

rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout

qsub -hold_jid "cs-${M}ethnicity-${ethnic1}_from_single" -l h=${skipnode} ALL.sh

