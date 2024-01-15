#!/usr/bin/env bash
#$ -N ukbb_geno_SAS
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=5000G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

#awk '{print $2}' /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref/ref_bim.txt > /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref/ref_bim.snp

ethnic='SAS'

mkdir /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/unrelatedness/${ethnic}/

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/all_data/${ethnic}/allchr \
--keep /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/${ethnic}.all_id \
--extract /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref/ref_bim.snp \
--rm-dup exclude-all \
--make-bed \
--maf 0.4 \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/unrelatedness/${ethnic}/allchr

#/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
#--bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/unrelatedness/${ethnic}/allchr \
#--genome \
#--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/unrelatedness/${ethnic}/allchr_genome

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \
--bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/unrelatedness/${ethnic}/allchr \
--rel-cutoff \
--keep-allele-order \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/unrelatedness/${ethnic}/allchr_unrel
