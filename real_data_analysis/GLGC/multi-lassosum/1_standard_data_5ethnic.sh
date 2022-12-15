#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-048|compute-053|compute-054|compute-066|compute-068|compute-070|compute-076|compute-113)'

M='5'

ethnic1='AFR_AMR_EAS_SAS_EUR'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum
runr(){
    R CMD BATCH --no-save --no-restore "$1" 1_write_R_standard_data.R 1_write_R_standard_data_${M}ethnicity_${ethnic1}.Rout
}
runr "--args ethnic=c('AFR','AMR','EAS','SAS','EUR') "

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum/${ethnic1}/1_standard_data

rm *.o*
rm *.e*
rm *.pe*
rm *.po*
rm *.Rout

qsub -l h=${skipnode} -hold_jid "mdata-*_HDL,tv_*_HDL" HDL.sh
qsub -l h=${skipnode} -hold_jid "mdata-*_LDL,tv_*_LDL" LDL.sh
qsub -l h=${skipnode} -hold_jid "mdata-*_logTG,tv_*_logTG" logTG.sh
qsub -l h=${skipnode} -hold_jid "mdata-*_nonHDL,tv_*_nonHDL" nonHDL.sh
qsub -l h=${skipnode} -hold_jid "mdata-*_TC,tv_*_TC" TC.sh


