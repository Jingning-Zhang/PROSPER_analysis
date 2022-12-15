

####################################################

library(dplyr)
library(readr)
library(stringr)

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data")
for(ethnic in c("AFR","AMR","EUR")){
  df_beta <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/height/tuning+validation/",ethnic,"_tuning.txt"))
  df_beta <- df_beta[,c(1,1,2)]
  write_tsv(df_beta, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_",ethnic,"_tuning.txt"), col_names=F)

  df_beta <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/height/tuning+validation/",ethnic,"_validation.txt"))
  df_beta <- df_beta[,c(1,1,2)]
  write_tsv(df_beta, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/pheno_",ethnic,"_testing.txt"), col_names=F)

  covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",ethnic,"_tuning.txt"))
  covar <- covar[,c(1,1:13)]
  colnames(covar) = c('fid','iid','sex','age',paste0('pc',1:10))
  write_tsv(covar, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_",ethnic,"_tuning.txt"))

  covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",ethnic,"_validation.txt"))
  covar <- covar[,c(1,1:13)]
  colnames(covar) = c('fid','iid','sex','age',paste0('pc',1:10))
  write_tsv(covar, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/sample_data/covar_",ethnic,"_testing.txt"))

}

