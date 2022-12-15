

####################################################

library(dplyr)
library(readr)
library(stringr)


EUR <- NULL
AFR <- NULL
AMR <- NULL
EAS <- NULL
SAS <- NULL
for(chr in 1:22){
  print(chr)

  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/EUR/rho_1_size_4_rep_1_GA_1_matched_to_ref_chr",chr,".RData"))
  df_beta$RSID <- unlist(lapply(str_split(df_beta$rsid,":"), FUN=function (x){x[1]}))
  tmp <- df_beta[,c("RSID","chr","a1","a0","beta","beta_se","n_eff")]
  colnames(tmp) <- c("rsid","chr","a1","a0","beta","beta_se","n_eff")
  EUR <- rbind(EUR, tmp)

  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/AFR/rho_1_size_1_rep_1_GA_1_matched_to_ref_chr",chr,".RData"))
  df_beta$RSID <- unlist(lapply(str_split(df_beta$rsid,":"), FUN=function (x){x[1]}))
  tmp <- df_beta[,c("RSID","chr","a1","a0","beta","beta_se","n_eff")]
  colnames(tmp) <- c("rsid","chr","a1","a0","beta","beta_se","n_eff")
  AFR <- rbind(AFR, tmp)

  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/AMR/rho_1_size_1_rep_1_GA_1_matched_to_ref_chr",chr,".RData"))
  df_beta$RSID <- unlist(lapply(str_split(df_beta$rsid,":"), FUN=function (x){x[1]}))
  tmp <- df_beta[,c("RSID","chr","a1","a0","beta","beta_se","n_eff")]
  colnames(tmp) <- c("rsid","chr","a1","a0","beta","beta_se","n_eff")
  AMR <- rbind(AMR, tmp)

  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/EAS/rho_1_size_1_rep_1_GA_1_matched_to_ref_chr",chr,".RData"))
  df_beta$RSID <- unlist(lapply(str_split(df_beta$rsid,":"), FUN=function (x){x[1]}))
  tmp <- df_beta[,c("RSID","chr","a1","a0","beta","beta_se","n_eff")]
  colnames(tmp) <- c("rsid","chr","a1","a0","beta","beta_se","n_eff")
  EAS <- rbind(EAS, tmp)

  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/SAS/rho_1_size_1_rep_1_GA_1_matched_to_ref_chr",chr,".RData"))
  df_beta$RSID <- unlist(lapply(str_split(df_beta$rsid,":"), FUN=function (x){x[1]}))
  tmp <- df_beta[,c("RSID","chr","a1","a0","beta","beta_se","n_eff")]
  colnames(tmp) <- c("rsid","chr","a1","a0","beta","beta_se","n_eff")
  SAS <- rbind(SAS, tmp)

}

write_tsv(EUR, "/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_summdata/EUR.txt")
write_tsv(AFR, "/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_summdata/AFR.txt")
write_tsv(AMR, "/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_summdata/AMR.txt")
write_tsv(EAS, "/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_summdata/EAS.txt")
write_tsv(SAS, "/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_summdata/SAS.txt")


#library(readr)
#for (l in 1:M){
#  if(ethnic[l]=="EUR"){
#    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_1_lassosum2/",ethnic[l],"/rep_1/best_tuning_rho_1_size_",4,".RData"))
#  }else{
#    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_1_lassosum2/",ethnic[l],"/rep_1/best_tuning_rho_1_size_1.RData"))
#  }
#  delta0 <- params2$delta[best]
#  lambda0 <- params2$lambda[best]
#  save(delta0,lambda0,file=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/best_lassosum_param/",ethnic[l],".RData"))
#}

