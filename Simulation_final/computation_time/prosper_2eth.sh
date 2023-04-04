#!/usr/bin/env bash
#$ -N prosper_2eth
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G,h='compute-119'
#$ -m e
#$ -t 1-10
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prosper/2pop/

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prosper/2pop/$SGE_TASK_ID

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge

date '+%A %W %Y %X'

Rscript lassosum2.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PRS-epr \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prosper/2pop/$SGE_TASK_ID/lassosum2 \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 22 \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EUR/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EUR/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam \
--testing F \
--NCORES 1

Rscript PRS-epr.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PRS-epr \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prosper/2pop/$SGE_TASK_ID/prosper \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 22 \
--lassosum_param /dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prosper/2pop/$SGE_TASK_ID/lassosum2/EUR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prosper/2pop/$SGE_TASK_ID/lassosum2/AFR/optimal_param.txt \
--NCORES 1

Rscript tuning_testing.R \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prosper/2pop/$SGE_TASK_ID/prosper \
--prefix AFR \
--testing F \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam \
--NCORES 1

date '+%A %W %Y %X'
