

######## new lassosum2
# SP method
rm(list=ls())

library(dplyr)

setting=1

mulmethods <- c(
                #"multi-lassosum",
                #"multi-lassosum(fixEUR)3",
                #"multi-lassosum_train_from_single_eth(fixEUR)",
                "multi-lassosum_train_from_single_eth"
                )

r2.vec <- numeric()
eth.vec <- character()
l_vec <- integer()
m_vec <- integer()
method_vec <- character()
ga_vec <- integer()
i <- 0
for (ethnic in c("AFR","AMR","EAS","SAS")){
  #ethnic='AFR'
for (size in 1:4){
  for (rho in 1:3){
      for (method in  c(mulmethods,
                        "lassosum2"

                        )){
        i <- i+1
        if(method %in% mulmethods){ mulflag <- T }else{ mulflag <- F }

        if(mulflag){

          if(method == "multi-lassosum_train_from_single_eth"){
            tryCatch({
              R2tmp <- tibble()
              for(rep in 1:10){
                tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/2_multi-lassosum_by_chr_10/Results/test_all_prs/",ethnic,"/rho_",rho,"_size_",size,"_validation_R2.rds"))
                R2tmp <- rbind(R2tmp, tmp)
              }
              R2 <- apply(as.matrix(R2tmp),MARGIN = 2, mean)
            }, error=function(e){})
          }
            a <- tryCatch({r2.vec[i:(i+5)] <- R2 }, error=function(e){})
            if(is.null(a)){r2.vec[i:(i+5)] <- NA}
            eth.vec[i:(i+5)] <- ethnic
            l_vec[i:(i+5)] <- rho
            m_vec[i:(i+5)] <- size
            method_vec[i:(i+5)] <- paste0(method," (", c("best_R2","best_R2_across_ancestry","lasso_combined_R2","lasso_combined_R2_cross_ancestry","superlearning_combined_R2","superlearning_combined_R2_cross_ancestry"),")")
            ga_vec[i:(i+5)] <- setting
            i <- i+5

          rm(R2)
        }

        if(!mulflag){
          if(method == "lassosum2"){
            tryCatch({
              R2tmp <- numeric()
              for(rep in 1:10){
                tmp <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,".rds"))
                R2tmp <- c(R2tmp, tmp)
              }
              R2 <- mean(R2tmp)
            }, error=function(e){})
          }

          a <- tryCatch({r2.vec[i] <- mean(R2) }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          l_vec[i] <- rho
          m_vec[i] <- size
          method_vec[i] <- method
          ga_vec[i] <- setting

          rm(R2)
        }

      }
    }
  }
}


res <- data.frame(eth.vec=eth.vec,l_vec=l_vec,m_vec=m_vec,r2.vec=r2.vec,method_vec=method_vec,ga_vec=ga_vec)
head(res)

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/result_table")
readr::write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/result_table/Restuls_GA_",setting,".txt"))
paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/result_table/Restuls_GA_",setting,".txt")


