#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge

date '+%A %W %Y %X'
#Monday 46 2022 01:19:08 AM

Rscript run_lassosum2.R \
--ref_dir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data \
--chrom 22 \
--sst_file /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AFR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AMR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EAS.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/SAS.txt \
--pop EUR,AFR,AMR,EAS,SAS \
--out_dir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2_chr22 \
--verbose 1

date '+%A %W %Y %X'
#Monday 46 2022 01:20:45 AM

Rscript tuning_testing.R \
--chrom 22 \
--pop EUR,AFR,AMR,EAS,SAS \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--out_dir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2_chr22 \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EUR/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AMR/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EAS/tuning_geno,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/SAS/tuning_geno \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EUR/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AMR/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EAS/pheno.fam,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/SAS/pheno.fam \
--testing FALSE \
--verbose 1
date '+%A %W %Y %X'
#Monday 46 2022 01:21:03 AM


Rscript run_PRS-epr.R \
--ref_dir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data \
--chrom 22 \
--sst_file /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EUR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AFR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/AMR.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/EAS.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/SAS.txt \
--pop EUR,AFR,AMR,EAS,SAS \
--lassosum_param /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2_chr22/EUR/best_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2_chr22/AFR/best_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2_chr22/AMR/best_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2_chr22/EAS/best_param.txt,/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/lassosum2_chr22/SAS/best_param.txt \
--out_dir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/PRSepr_chr22 \
--verbose 1

date '+%A %W %Y %X'
#Monday 46 2022 01:22:58 AM

Rscript ensemble_tuning_testing.R \
--chrom 22 \
--PATH_plink /dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \
--out_dir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/5pop/PRSepr_chr22 \
--testing TRUE \
--bfile_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno \
--pheno_tuning /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam \
--bfile_testing /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/testing_geno \
--pheno_testing /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam \
--verbose 1

date '+%A %W %Y %X'
#Monday 46 2022 01:23:44 AM


