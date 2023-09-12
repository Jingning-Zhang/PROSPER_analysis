# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

setting=1
rep=1

para='10'

M='5'

ethnic1="AFR_AMR_EAS_EUR_SAS"
ethnic=c("AFR","AMR","EAS","EUR","SAS")

L=para
Lc=para

rho_vec <- integer()
size_vec <- integer()
eth_vec <- character()
Ntuning_vec <- integer()
R2 <- numeric()
for (rho in 1:3){
  for (size in 1:4){
    for (targetethnic in ethnic){
      for (Ntuning in c(5000,3000,1000,500,300,100)){
        res <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/number_of_tuning_samples/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para,"/Results/test_all_prs/",targetethnic,"/rho_",rho,"_size_",size,"_Ntuning",Ntuning,"_validation_R2.rds"))[[1]]
        rho_vec <- c(rho_vec, rho)
        size_vec <- c(size_vec, size)
        eth_vec <- c(eth_vec, targetethnic)
        Ntuning_vec <- c(Ntuning_vec, Ntuning)
        R2 <- c(R2, res)
      }
    }
  }
}
res <- data.frame(rho_vec, size_vec, eth_vec, Ntuning_vec, R2)
library(readr)
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/number_of_tuning_samples/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para,"/R2.txt"))
