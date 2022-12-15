#!/usr/bin/env bash
#$ -N pip_AFR
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/
mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/ref_geno/
for ethnic in "EUR" "AFR" "AMR"
do
  mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/ref_geno/${ethnic}
for trait in "bmi" "height"
do
  mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/ref_geno/${ethnic}/${trait}/
for chr in {1..22}
do
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--threads 1 \
--bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh38/${ethnic}/chr${chr} \
--extract /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/snps/${ethnic}_${trait}_rsid.txt \
--force-intersect \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/ref_geno/${ethnic}/${trait}/ref_chr${chr}
done
done
done

#/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
#--threads 1 \
#--bfile /dcs04/nilanjan/data/jzhang2/UKBB/genotype/tuning+validation/tuning/EUR/chr1 \
#--freq \
#--out /users/jzhang2/tmp_tmp/test
