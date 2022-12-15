#!/usr/bin/env bash
#$ -N prsepr_2eth
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=22G,h='(compute-119|compute-124|compute-129)'
#$ -m e
#$ -t 1-10
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge

date '+%A %W %Y %X'

Rscript lassosum2.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2 \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/AFR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/AMR.txt \
--pop EUR,AFR,AMR \
--chrom 1-22 \
--Ll 10 --Ld 10 \
--testing T \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/EUR/allchr,/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/AFR/allchr,/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/AMR/allchr \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_EUR_tuning.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_AFR_tuning.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_AMR_tuning.txt \
--covar_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_EUR_tuning.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_AFR_tuning.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_AMR_tuning.txt \
--bfile_testing /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/EUR/allchr,/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/AFR/allchr,/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/AMR/allchr \
--pheno_testing /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_EUR_testing.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_AFR_testing.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_AMR_testing.txt \
--covar_testing /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_EUR_testing.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_AFR_testing.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_AMR_testing.txt \
--NCORES 22


Rscript PRS-epr.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/PRSepr \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/AFR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/AMR.txt \
--pop EUR,AFR,AMR \
--chrom 1-22 \
--Ll 10 --Lc 10 \
--lassosum_param /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2/EUR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2/AFR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2/AMR/optimal_param.txt \
--NCORES 22

Rscript tuning_testing.R \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/PRSepr \
--prefix AFR \
--testing T \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/AFR/allchr \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_AFR_tuning.txt \
--covar_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_AFR_tuning.txt \
--bfile_testing /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/AFR/allchr \
--pheno_testing /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_AFR_testing.txt \
--covar_testing /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_AFR_testing.txt \
--NCORES 100



## parameters from original study

Rscript PRS-epr.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/PRSepr_from_original_study \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/AFR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/AMR.txt \
--pop EUR,AFR,AMR \
--chrom 1-22 \
--Ll 10 --Lc 10 \
--lassosum_param /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2_from_original_study/EUR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2_from_original_study/AFR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2_from_original_study/AMR/optimal_param.txt \
--NCORES 22

Rscript tuning_testing.R \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/PRSepr_from_original_study \
--prefix AFR \
--testing T \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/AFR/allchr \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_AFR_tuning.txt \
--covar_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_AFR_tuning.txt \
--bfile_testing /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/AFR/allchr \
--pheno_testing /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_AFR_testing.txt \
--covar_testing /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_AFR_testing.txt \
--NCORES 100
