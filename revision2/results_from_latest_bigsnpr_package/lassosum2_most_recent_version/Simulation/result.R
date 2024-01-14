
rm(list=ls())

library(readr)

rep=3

setting_vec <- integer()
rho_vec <- integer()
size_vec <- integer()
eth_vec <- character()
PROSPER <- numeric()
lassosum2 <- numeric()
lassosum2_EUR <- numeric()
lassosum2_sl <- numeric()
lassosum2_wt <- numeric()
lassosum2_old <- numeric()
lassosum2_wt_old <- numeric()
for(setting in 1:5){
  for(rho in 1:3){
    for(size in 1:4){
      for(eth in c('AFR','AMR','EAS','SAS')){

        a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/2_multi-lassosum_by_chr_10/Results/test_all_prs/",eth,"/rho_",rho,"_size_",size,"_validation_R2.rds"))$superlearning_combined_R2_cross_ancestry  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; PROSPER <- c(PROSPER,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",eth,"/rep_",rep,"/lassosum2/prs/results/rho_",rho,"_size_",size,"/",eth,"_validation_R2_best.rds"))$best_R2  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; lassosum2 <- c(lassosum2,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",eth,"/rep_",rep,"/lassosum2/prs/results/rho_",rho,"_size_",size,"/",eth,"_validation_R2_EUR_best.rds"))$EUR_R2  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; lassosum2_EUR <- c(lassosum2_EUR,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",eth,"/rep_",rep,"/lassosum2/prs/results/rho_",rho,"_size_",size,"/",eth,"_validation_R2_sl.rds"))$superlearning_combined_R2  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; lassosum2_sl <- c(lassosum2_sl,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",eth,"/rep_",rep,"/lassosum2/prs/results/rho_",rho,"_size_",size,"/",eth,"_validation_R2_sl.rds"))$weighted_combined_R2 }, error=function(e){})
        if(is.null(a)){tmp <- NA}; lassosum2_wt <- c(lassosum2_wt,tmp)

        #a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",eth,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,".rds")) }, error=function(e){})
        #if(is.null(a)){tmp <- NA}; lassosum2_old <- c(lassosum2_old,tmp)

        #a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",eth,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,"_weighted_lassosum2.rds"))$R2_5eth }, error=function(e){})
        #if(is.null(a)){tmp <- NA}; lassosum2_wt_old <- c(lassosum2_wt_old,tmp)

        setting_vec <- c(setting_vec, setting)
        rho_vec <- c(rho_vec, rho)
        size_vec <- c(size_vec, size)
        eth_vec <- c(eth_vec, eth)
      }
    }
  }
}


res <- data.frame(setting=setting_vec,
                  rho=rho_vec,
                  size=size_vec,
                  ethnic=eth_vec,
                  lassosum2=lassosum2,
                  lassosum2_EUR=lassosum2_EUR,
                  #lassosum2_old = lassosum2_old,
                  lassosum2_wt = lassosum2_wt,
                  #lassosum2_wt_old = lassosum2_wt_old,
                  lassosum2_sl=lassosum2_sl,
                  PROSPER=PROSPER
)
sum(is.na(res))

write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Simulation_lassosum2_sl_rev2_rep",rep,".txt"))


for(rep in 1:3){
  tmp <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Simulation_lassosum2_sl_rev2_rep",rep,".txt"))
  tmp0 <- tmp[,1:4]
  tmp1 <- tmp[,-1:-4]
  if(rep==1){res <- tmp1}else{res <- res+tmp1}
}
res <- cbind(tmp0, res/3)
write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Simulation_lassosum2_sl_rev2_rep_ave.txt"))
