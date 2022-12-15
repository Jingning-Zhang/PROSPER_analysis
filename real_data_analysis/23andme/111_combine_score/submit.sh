#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-048|compute-053|compute-054|compute-066|compute-068|compute-070|compute-076|compute-113)'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/23andme/111_combine_score/9_compress_score_+weighted
rm *.o*
rm *.e*
rm *.pe*
rm *.po*
rm *.Rout

qsub -l h=${skipnode} -hold_jid "cs23-*" ALL.sh



