
rm(list=ls())

library(readr)

trait_vec <- integer()
eth_vec <- character()
PROSPER <- numeric()
lassosum2 <- numeric()
lassosum2_EUR <- numeric()
lassosum2_sl <- numeric()
lassosum2_wt <- numeric()
lassosum2_old <- numeric()
lassosum2_wt_old <- numeric()
for(trait in c("height","bmi")){
      for(eth in c('AFR')){

        a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum_train_from_single_eth/AFR_AMR_EUR/2_multi-lassosum_by_chr_10/Results/final_results/",eth,"/",trait,"_validation_R2.rds"))$superlearning_combined_R2_cross_ancestry  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; PROSPER <- c(PROSPER,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait,"/lassosum2/prs/results/",eth,"_validation_R2_best.rds"))$best_R2  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; lassosum2 <- c(lassosum2,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait,"/lassosum2/prs/results/",eth,"_validation_R2_EUR_best.rds"))$EUR_R2  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; lassosum2_EUR <- c(lassosum2_EUR,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait,"/lassosum2/prs/results/",eth,"_validation_R2_sl.rds"))$superlearning_combined_R2  }, error=function(e){})
        if(is.null(a)){tmp <- NA}; lassosum2_sl <- c(lassosum2_sl,tmp)

        a <- tryCatch({ tmp <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait,"/lassosum2/prs/results/",eth,"_validation_R2_sl.rds"))$weighted_combined_R2 }, error=function(e){})
        if(is.null(a)){tmp <- NA}; lassosum2_wt <- c(lassosum2_wt,tmp)

        #a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",eth,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,".rds")) }, error=function(e){})
        #if(is.null(a)){tmp <- NA}; lassosum2_old <- c(lassosum2_old,tmp)

        #a <- tryCatch({ tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",eth,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,"_weighted_lassosum2.rds"))$R2_5eth }, error=function(e){})
        #if(is.null(a)){tmp <- NA}; lassosum2_wt_old <- c(lassosum2_wt_old,tmp)

        trait_vec <- c(trait_vec, trait)
        eth_vec <- c(eth_vec, eth)

  }
}


res <- data.frame(trait=trait_vec,
                  ethnic=eth_vec,
                  lassosum2=lassosum2,
                  lassosum2_EUR=lassosum2_EUR,
                  #lassosum2_old = lassosum2_old,
                  lassosum2_wt = lassosum2_wt,
                  #lassosum2_wt_old = lassosum2_wt_old,
                  lassosum2_sl=lassosum2_sl,
                  PROSPER=PROSPER
)

mean(res$PROSPER/res$lassosum2_sl)
#1.117314

write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/AoU_lassosum2_sl_rev2.txt"))

