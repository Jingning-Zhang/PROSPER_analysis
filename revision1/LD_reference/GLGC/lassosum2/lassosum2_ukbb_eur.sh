#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-048|compute-053|compute-054|compute-057|compute-066|compute-068|compute-070|compute-076|compute-113)'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/lassosum2/1_prepare_data_lassosum2
qsub -l h=${skipnode} -hold_jid ukbb_ld_eur_merge EUR_HDL.sh

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/lassosum2/2_run_lassosum2
rm *.e*
rm *.o*
rm *.Rout
qsub -l h=${skipnode} -hold_jid data_EUR_HDL EUR_HDL.sh

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/lassosum2/3_lassosum2_clean_score
rm *.e*
rm *.o*
rm *.Rout
qsub -l h=${skipnode} -hold_jid run_EUR_HDL EUR_HDL.sh

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/lassosum2/4_ukbb_tuning+validation
rm *.e*
rm *.o*
rm *.Rout
qsub -l h=${skipnode} -hold_jid run_EUR_HDL EUR_HDL.sh


