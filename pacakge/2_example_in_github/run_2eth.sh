#!/usr/bin/env bash
#$ -N prsepr_2eth
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=100G,h='(compute-119|compute-124|compute-129)'
#$ -m e
#$ -t 1-10
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

package='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PRS-epr'
path_example='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/example_in_github'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge

date '+%A %W %Y %X'

Rscript lassosum2.R \
--PATH_package ${package} \
--PATH_out ${path_example}/PRS-epr_example_results/lassosum2 \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_example}/sample_data/EUR/tuning_geno,${path_example}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--bfile_testing ${path_example}/sample_data/EUR/testing_geno,${path_example}/sample_data/AFR/testing_geno \
--pheno_testing ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--testing T \
--NCORES 22

Rscript PRS-epr.R \
--PATH_package ${package} \
--PATH_out ${path_example}/PRS-epr_example_results/PRSepr \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--lassosum_param ${path_example}/PRS-epr_example_results/lassosum2/EUR/optimal_param.txt,${path_example}/PRS-epr_example_results/lassosum2/AFR/optimal_param.txt \
--chrom 1-22 \
--NCORES 22

Rscript tuning_testing.R \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--PATH_out ${path_example}/PRS-epr_example_results/PRSepr \
--prefix AFR \
--testing T \
--bfile_tuning ${path_example}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_example}/sample_data/AFR/pheno.fam \
--bfile_testing ${path_example}/sample_data/AFR/testing_geno \
--pheno_testing ${path_example}/sample_data/AFR/pheno.fam \
--NCORES 22

date '+%A %W %Y %X'


