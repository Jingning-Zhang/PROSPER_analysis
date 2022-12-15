

####################################################

library(dplyr)
library(readr)
library(stringr)

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata")
for(ethnic in c("AFR","AMR","EUR")){
  df_beta <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/",ethnic,"_height_df_beta.rds"))
  df_beta <- df_beta[,c("rsid","chr","a1","a0","beta","beta_se","n_eff")]
  write_tsv(df_beta, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/summdata/",ethnic,".txt"))

}
