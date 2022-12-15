#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

for ethnic in "AFR" "AMR" "EAS" "EUR" "SAS"
do

cd /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data/${ethnic}

for (( i = 1; i < 23; i++ )); do
cp /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block/${ethnic}/standard_data/chr${i}_LD.RData ./
done

done

#/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir/plink2 \
#  --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block/AFR/tmp/byblock/chr1/chr1_247344518_249239466 \
#  --freq \
#  --out /users/jzhang2/tmp_tmp/freq