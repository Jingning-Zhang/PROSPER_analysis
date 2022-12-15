#!/usr/bin/env bash
#$ -N prsepr_2eth
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -m e
#$ -pe local 5
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/2pop

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge

date '+%A %W %Y %X'

Rscript lassosum2.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PRS-epr \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/2pop/lassosum2 \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EUR/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EUR/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam \
--testing F \
--NCORES 5

Rscript PRS-epr.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PRS-epr \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/2pop/PRSepr \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--lassosum_param /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/2pop/lassosum2/EUR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/2pop/lassosum2/AFR/optimal_param.txt \
--NCORES 5

Rscript tuning_testing.R \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/2pop/PRSepr \
--prefix AFR \
--testing F \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam \
--NCORES 5

date '+%A %W %Y %X'
