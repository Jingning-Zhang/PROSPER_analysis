#!/usr/bin/env bash
#$ -N LDref
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu


mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/
mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno

for ethnic in "EUR" "AFR" "AMR" "EAS" "SAS"
do
  mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/${ethnic}

for chr in {1..22}
do
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--threads 1 \
--bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/${ethnic}/chr${chr} \
--extract /dcl01/chatterj/data/jin/MEGA/mega-hm3-rsid.txt \
--rm-dup exclude-all \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/${ethnic}/ref_chr${chr}
done
done


for ethnic in "EUR" "AFR" "AMR" "EAS" "SAS"
do
for chr in {1..22}
do
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--threads 1 \
--bfile /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/${ethnic}/ref_chr${chr} \
--freq \
--out /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/${ethnic}/freq_ref_chr${chr}
done
done

