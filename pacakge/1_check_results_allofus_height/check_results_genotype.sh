#!/usr/bin/env bash
#$ -N prsepr_2eth
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G,h='(compute-119|compute-124|compute-129)'
#$ -m e
#$ -t 1-10
#$ -M jzhan218@jhu.edu

for ethnic in "AFR" "AMR" "EAS" "EUR" "SAS"
do

for i in {1..22}
do
echo /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/${ethnic}/chr${i} >> /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/${ethnic}/allchr_mergelist
done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--keep-allele-order \
--make-bed \
--merge-list /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/${ethnic}/allchr_mergelist \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/${ethnic}/allchr \
--threads 100

done


for ethnic in "AFR" "AMR" "EAS" "EUR" "SAS"
do

for i in {1..22}
do
echo /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/${ethnic}/chr${i} >> /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/${ethnic}/allchr_mergelist
done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--keep-allele-order \
--make-bed \
--merge-list /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/${ethnic}/allchr_mergelist \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/${ethnic}/allchr \
--threads 100

done
