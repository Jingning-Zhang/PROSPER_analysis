#!/usr/bin/env bash
#$ -N LDref
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu


mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/
for ethnic in "EUR" "AFR" "AMR" "EAS" "SAS"
do
  mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/${ethnic}
for trait in "HDL" "LDL" "logTG" "nonHDL" "TC"
do
  mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/${ethnic}/${trait}/
for chr in {1..22}
do
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--threads 1 \
--bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/${ethnic}/chr${chr} \
--extract /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/snps/${ethnic}_${trait}_rsid.txt \
--force-intersect \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/${ethnic}/${trait}/ref_chr${chr}
done
done
done
