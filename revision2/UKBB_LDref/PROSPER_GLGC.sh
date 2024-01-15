#!/usr/bin/env bash
#$ -N GLGC
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -m e
#$ -t 1-4
#$ -M jzhan218@jhu.edu

#module load conda_R/4.0

package='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref'
path_example='/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/compare_lassosum_and_prosper/GLGC'
ukbb_geno_path='/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype'


cd /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/try_from_github/PROSPER/scripts

traits=('HDL' 'LDL' 'logTG' 'TC')
#trait=${traits[$SGE_TASK_ID]}

trait=${traits[1]}

mkdir ${path_example}/
mkdir ${path_example}/results_ukbref/
mkdir ${path_example}/results_ukbref/${trait}
mkdir ${path_example}/results_ukbref/${trait}/lassosum2
mkdir ${path_example}/results_ukbref/${trait}/PROSPER

date '+%A %W %Y %X'

Rscript lassosum2.R \
--PATH_package ${package} \
--PATH_out ${path_example}/results_ukbref/${trait}/lassosum2 \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--FILE_sst ${path_example}/summdata/EUR_${trait}.txt,${path_example}/summdata/AFR_${trait}.txt,${path_example}/summdata/AMR_${trait}.txt,${path_example}/summdata/EAS_${trait}.txt,${path_example}/summdata/SAS_${trait}.txt \
--pop EUR,AFR,AMR,EAS,SAS \
--chrom 1-22 \
--bfile_tuning ${ukbb_geno_path}/tuning/EUR/allchr,${ukbb_geno_path}/tuning/AFR/allchr,${ukbb_geno_path}/tuning/AMR/allchr,${ukbb_geno_path}/tuning/EAS/allchr,${ukbb_geno_path}/tuning/SAS/allchr \
--pheno_tuning ${path_example}/pheno_data/pheno_EUR_${trait}_tuning.txt,${path_example}/pheno_data/pheno_AFR_${trait}_tuning.txt,${path_example}/pheno_data/pheno_AMR_${trait}_tuning.txt,${path_example}/pheno_data/pheno_EAS_${trait}_tuning.txt,${path_example}/pheno_data/pheno_SAS_${trait}_tuning.txt \
--covar_tuning ${path_example}/pheno_data/covar_EUR_${trait}_tuning.txt,${path_example}/pheno_data/covar_AFR_${trait}_tuning.txt,${path_example}/pheno_data/covar_AMR_${trait}_tuning.txt,${path_example}/pheno_data/covar_EAS_${trait}_tuning.txt,${path_example}/pheno_data/covar_SAS_${trait}_tuning.txt \
--bfile_testing ${ukbb_geno_path}/validation/EUR/allchr,${ukbb_geno_path}/validation/AFR/allchr,${ukbb_geno_path}/validation/AMR/allchr,${ukbb_geno_path}/validation/EAS/allchr,${ukbb_geno_path}/validation/SAS/allchr \
--pheno_testing ${path_example}/pheno_data/pheno_EUR_${trait}_testing.txt,${path_example}/pheno_data/pheno_AFR_${trait}_testing.txt,${path_example}/pheno_data/pheno_AMR_${trait}_testing.txt,${path_example}/pheno_data/pheno_EAS_${trait}_testing.txt,${path_example}/pheno_data/pheno_SAS_${trait}_testing.txt \
--covar_testing ${path_example}/pheno_data/covar_EUR_${trait}_testing.txt,${path_example}/pheno_data/covar_AFR_${trait}_testing.txt,${path_example}/pheno_data/covar_AMR_${trait}_testing.txt,${path_example}/pheno_data/covar_EAS_${trait}_testing.txt,${path_example}/pheno_data/covar_SAS_${trait}_testing.txt \
--testing T \
--cleanup F \
--NCORES 22

Rscript PROSPER.R \
--PATH_package ${package} \
--PATH_out ${path_example}/results_ukbref/${trait}/PROSPER \
--FILE_sst ${path_example}/summdata/EUR_${trait}.txt,${path_example}/summdata/AFR_${trait}.txt,${path_example}/summdata/AMR_${trait}.txt,${path_example}/summdata/EAS_${trait}.txt,${path_example}/summdata/SAS_${trait}.txt \
--pop EUR,AFR,AMR,EAS,SAS \
--lassosum_param ${path_example}/results_ukbref/${trait}/lassosum2/EUR/optimal_param.txt,${path_example}/results_ukbref/${trait}/lassosum2/AFR/optimal_param.txt,${path_example}/results_ukbref/${trait}/lassosum2/AMR/optimal_param.txt,${path_example}/results_ukbref/${trait}/lassosum2/EAS/optimal_param.txt,${path_example}/results_ukbref/${trait}/lassosum2/SAS/optimal_param.txt \
--chrom 1-22 \
--NCORES 22

for eth in "EUR" "AFR" "AMR" "EAS" "SAS"
do

Rscript tuning_testing.R \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--PATH_out ${path_example}/results_ukbref/${trait}/PROSPER \
--prefix ${eth} \
--testing T \
--cleanup F \
--bfile_tuning ${ukbb_geno_path}/tuning/${eth}/allchr \
--pheno_tuning ${path_example}/pheno_data/pheno_${eth}_${trait}_tuning.txt \
--covar_tuning ${path_example}/pheno_data/covar_${eth}_${trait}_tuning.txt \
--bfile_testing ${ukbb_geno_path}/validation/${eth}/allchr \
--pheno_testing ${path_example}/pheno_data/pheno_${eth}_${trait}_testing.txt \
--covar_testing ${path_example}/pheno_data/covar_${eth}_${trait}_testing.txt \
--NCORES 22

done

date '+%A %W %Y %X'
