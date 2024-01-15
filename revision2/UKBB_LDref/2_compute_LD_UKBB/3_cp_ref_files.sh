#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref

for ethnic in "AFR" "EAS" "SAS"
do

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref/${ethnic}

cd /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref/${ethnic}

for (( i = 1; i < 23; i++ )); do
cp /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block_UKBB/${ethnic}/standard_data/chr${i}_LD.RData ./
done

done

