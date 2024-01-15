#!/usr/bin/env bash
#$ -N ukbb_ld_EUR
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/
mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/EUR/

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/EUR/tmp
awk '{print $2}' /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/EUR/chr$SGE_TASK_ID.bim > /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/EUR/tmp/1000g_chr$SGE_TASK_ID.snp

/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir/plink2 \
  --pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr$SGE_TASK_ID \
  --extract-intersect /dcl01/chatterj/data/jin/MEGA/mega-hm3-rsid.txt /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/EUR/tmp/1000g_chr$SGE_TASK_ID.snp \
  --keep /dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/unrelated_whites.id \
  --ref-allele /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/EUR/chr$SGE_TASK_ID.bim 6 2 \
  --rm-dup exclude-mismatch \
  --make-bed \
  --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/EUR/ref_chr$SGE_TASK_ID \
  --threads 1

