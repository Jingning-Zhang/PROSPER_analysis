#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

M='3'

ethnic1='AFR_AMR_EUR'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/multi-lassosum
runr(){
    R CMD BATCH --no-save --no-restore "$1" 1_write_R_standard_data.R 1_write_R_standard_data_${M}ethnicity_${ethnic1}.Rout
}
runr "--args ethnic=c('AFR','AMR','EUR') "

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/multi-lassosum/${ethnic1}/1_standard_data
rm *.o*
rm *.e*
rm *.pe*
rm *.po*
rm *.Rout

qsub -l h=${skipnode} -hold_jid "mdata-*_height,tv_*_height" height.sh
qsub -l h=${skipnode} -hold_jid "mdata-*_bmi,tv_*_bmi" bmi.sh


