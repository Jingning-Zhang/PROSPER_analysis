#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-048|compute-053|compute-054|compute-066|compute-068|compute-070|compute-076|compute-113)'


para='10'
M='5'

ethnic1='AFR_AMR_EAS_SAS_EUR'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth
runr(){
    R CMD BATCH --no-save --no-restore "$1" 3_write_R_multi-lassosum_clean_score.R 3_write_R_multi-lassosum_clean_score_${M}ethnicity_${ethnic1}.Rout
}
runr "--args ethnic=c('AFR','AMR','EAS','SAS','EUR') L='${para}' Lc='${para}' para='${para}'"


cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/${ethnic1}/3_multi-lassosum_by_chr_${para}_clean_score

rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout

qsub -l h=${skipnode} -hold_jid "meprs-${M}ethnicity-${ethnic1}" ALL.sh




