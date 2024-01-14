
rm(list=ls())

library(readr)

trait_vec <- integer()
eth_vec <- character()
PROSPER <- numeric()
ldpred2 <- numeric()
ldpred2_EUR <- numeric()
ldpred2_sl <- numeric()
ldpred2_wt <- numeric()
ldpred2_old <- numeric()
ldpred2_wt_old <- numeric()
for(trait in c("height","bmi")){
      for(eth in c('AFR')){

        a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum_train_from_single_eth/AFR_AMR_EUR/2_multi-lassosum_by_chr_10/Results/final_results/",eth,"/",trait,"_validation_R2.rds"))$superlearning_combined_R2_cross_ancestry  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; PROSPER <- c(PROSPER,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait,"/ldpred2/prs/results/",eth,"_validation_R2_best.rds"))$best_R2  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; ldpred2 <- c(ldpred2,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait,"/ldpred2/prs/results/",eth,"_validation_R2_EUR_best.rds"))$best_R2  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; ldpred2_EUR <- c(ldpred2_EUR,tmp)

        #a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait,"/ldpred2/prs/results/",eth,"_validation_R2_sl.rds"))$superlearning_combined_R2  }, error=function(e){})
        #if(is.null(a)){tmp <- NA}; ldpred2_sl <- c(ldpred2_sl,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait,"/ldpred2/prs/results/",eth,"_allancs.rds"))$weighted_all }, error=function(e){})
        if(is.null(a)){tmp <- NA}; ldpred2_wt <- c(ldpred2_wt,tmp)

        #a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_ldpred2/",eth,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,".rds")) }, error=function(e){})
        #if(is.null(a)){tmp <- NA}; ldpred2_old <- c(ldpred2_old,tmp)

        #a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_ldpred2/",eth,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,"_weighted_ldpred2.rds"))$R2_5eth }, error=function(e){})
        #if(is.null(a)){tmp <- NA}; ldpred2_wt_old <- c(ldpred2_wt_old,tmp)

        trait_vec <- c(trait_vec, trait)
        eth_vec <- c(eth_vec, eth)

  }
}


res <- data.frame(trait=trait_vec,
                  ethnic=eth_vec,
                  ldpred2=ldpred2,
                  ldpred2_EUR=ldpred2_EUR,
                  #ldpred2_old = ldpred2_old,
                  ldpred2_wt = ldpred2_wt,
                  #ldpred2_wt_old = ldpred2_wt_old,
                  #ldpred2_sl=ldpred2_sl,
                  PROSPER=PROSPER
)

write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/AoU_ldpred2_rev2.txt"))


