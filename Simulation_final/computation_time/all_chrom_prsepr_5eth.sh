#!/usr/bin/env bash
#$ -N less
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -m e
#$ -pe local 5
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge

date '+%A %W %Y %X'

Rscript lassosum2.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PRS-epr \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2 \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AFR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AMR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EAS.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/SAS.txt \
--pop EUR,AFR,AMR,EAS,SAS \
--chrom 1-22 \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EUR/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AMR/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EAS/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/SAS/tuning_geno \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EUR/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AMR/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EAS/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/SAS/pheno.fam \
--testing F \
--NCORES 5

Rscript PRS-epr.R \
--PATH_package /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PRS-epr \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/PRSepr \
--FILE_sst /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AFR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AMR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EAS.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/SAS.txt \
--pop EUR,AFR,AMR,EAS,SAS \
--chrom 1-22 \
--lassosum_param /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2/EUR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2/AFR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2/AMR/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2/EAS/optimal_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2/SAS/optimal_param.txt \
--NCORES 5

Rscript tuning_testing.R \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--PATH_out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/PRSepr \
--prefix AFR \
--testing F \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam \
--NCORES 5

date '+%A %W %Y %X'

# compute-141

# -l mem_free=3G,h_vmem=3G,h_fsize=100G,h='(compute-119|compute-124|compute-129)'
