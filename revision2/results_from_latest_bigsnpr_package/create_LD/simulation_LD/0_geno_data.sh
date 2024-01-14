#!/usr/bin/env bash
#$ -N Sim
#$ -cwd
#$ -l mem_free=1G,h_vmem=1G,h_fsize=100G
#$ -m e
#$ -t 1-110
#$ -M jzhan218@jhu.edu

# run this inside ./0_geno_data

readarray -t races < ./race.para
readarray -t chrs < ./chr.para

race=${races[$(($SGE_TASK_ID-1))]}
chr=${chrs[$(($SGE_TASK_ID-1))]}

mkdir  /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/compare_lassosum_and_prosper/Simulation/geno_data/
mkdir  /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/compare_lassosum_and_prosper/Simulation/geno_data/${race}

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/${race}/tv_chr${chr} \
--keep /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/${race}/test.id.txt \
--update-name /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/compare_lassosum_and_prosper/Simulation/geno_data/rsid_map \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/compare_lassosum_and_prosper/Simulation/geno_data/${race}/tuning_chr${chr} \
--threads 1

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/${race}/tv_chr${chr} \
--keep /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/${race}/validation.id.txt \
--update-name /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/compare_lassosum_and_prosper/Simulation/geno_data/rsid_map \
--make-bed \
--out /dcs04/nilanjan/data/jzhang2/MEPRS/revision2/compare_lassosum_and_prosper/Simulation/geno_data/${race}/validation_chr${chr} \
--threads 1

