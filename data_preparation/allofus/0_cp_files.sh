#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -t 1-5
#$ -m e
#$ -M jzhan218@jhu.edu


cp -r /dcs04/nilanjan/data/rzhao/AllofUs/GWAS/plink_height /dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/allofus_raw/height
cp -r /dcs04/nilanjan/data/rzhao/AllofUs/GWAS/plink_BMI /dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/allofus_raw/bmi

cd /dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/allofus_raw/height


