#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-048|compute-053|compute-054|compute-066|compute-068|compute-076|compute-113)'


para='10'
M='3'

ethnic1='AFR_AMR_EUR'


cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/multi-lassosum_train_from_single_eth
runr(){
    R CMD BATCH --no-save --no-restore "$1" 2_write_R_multi-lassosum.R 2_write_R_multi-lassosum_${M}ethnicity_${ethnic1}.Rout
}
runr "--args ethnic=c('AFR','AMR','EUR') L='${para}' Lc='${para}' para='${para}'"

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/multi-lassosum_train_from_single_eth/${ethnic1}/2_multi-lassosum_by_chr_${para}

rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout

qsub -l h=${skipnode} -hold_jid "mdata-${M}ethnicity-${ethnic1}_*" ALL.sh

