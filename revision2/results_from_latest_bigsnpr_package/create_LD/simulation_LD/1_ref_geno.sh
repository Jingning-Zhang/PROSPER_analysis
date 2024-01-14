#!/usr/bin/env bash
#$ -N LDref
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu


for ethnic in "EUR" "AFR" "AMR" "EAS" "SAS"
do
for chr in {1..22}
do
/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--threads 1 \
--bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/${ethnic}/ref_chr${chr} \
--freq \
--out /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/${ethnic}/freq_ref_chr${chr}
done
done

