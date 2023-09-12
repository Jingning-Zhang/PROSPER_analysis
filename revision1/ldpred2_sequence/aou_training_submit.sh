#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-048|compute-053|compute-054|compute-066|compute-068|compute-070|compute-076|compute-113)'

ind1='1'
ind2='1'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth
runr(){
    R CMD BATCH --no-save --no-restore "$1" aou_training.R ${ind1}${ind2}-aou_training.Rout
}
runr "--args temp1='${ind1}' temp2='${ind2}'"

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/${ethnic1}/2_multi-lassosum_by_chr_${para}

rm *.e*
rm *.o*
rm *.pe*
rm *.po*
rm *.Rout

qsub -l h=${skipnode} -hold_jid "mdata-${M}ethnicity-${ethnic1}_*" ALL.sh

