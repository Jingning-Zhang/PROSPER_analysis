#!/bin/bash
#SBATCH --job-name ukbb_ld_AFR
#SBATCH --mem=10G
#SBATCH --array=1-22
#SBATCH --exclude=compute-070,compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/
mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/AFR/

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/AFR/tmp
awk '{print $2}' /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/AFR/chr$SLURM_ARRAY_TASK_ID.bim > /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/AFR/tmp/1000g_chr$SLURM_ARRAY_TASK_ID.snp

/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir/plink2 \
  --pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr$SLURM_ARRAY_TASK_ID \
  --extract-intersect /dcl01/chatterj/data/jin/MEGA/mega-hm3-rsid.txt /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/AFR/tmp/1000g_chr$SLURM_ARRAY_TASK_ID.snp \
  --keep /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/unrelatedness/AFR/allchr_unrel.rel.id \
  --ref-allele /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/AFR/chr$SLURM_ARRAY_TASK_ID.bim 6 2 \
  --rm-dup exclude-mismatch \
  --make-bed \
  --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/AFR/ref_chr$SLURM_ARRAY_TASK_ID \
  --threads 1

