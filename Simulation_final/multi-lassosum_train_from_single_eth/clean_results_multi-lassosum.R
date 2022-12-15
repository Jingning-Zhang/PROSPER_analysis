
######## lassosum2

rm(list=ls())

setting=2
rep=1

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
      for (method in  c("2-ethnic lassosum",
                        "2-ethnic lassosum (fix_EUR_coef)",
                        "5-ethnic lassosum",
                        "5-ethnic lassosum (fix_EUR_coef)",

                        "combined 2-ethnic lassosum",
                        "combined 2-ethnic lassosum (fix_EUR_coef)",
                        "combined 5-ethnic lassosum",
                        "combined 5-ethnic lassosum (fix_EUR_coef)",

                        "less_para combined 2-ethnic lassosum",
                        "less_para combined 2-ethnic lassosum (fix_EUR_coef)",
                        "less_para combined 5-ethnic lassosum",
                        "less_para combined 5-ethnic lassosum (fix_EUR_coef)",

                        "lassosum2",
                        "EUR lassosum2",
                        "Weighted lassosum2"
                        )){
        i <- i+1
        if(method == "lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,".rds")) }, error=function(e){})
        }else if(method == "EUR lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_4_EUR.rds")) }, error=function(e){})
        }else if(method == "Weighted lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,"_weighted_lassosum2.rds")) }, error=function(e){})
        }else if(method == "2-ethnic lassosum"){
          tryCatch({
            R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/2_multi-lassosum_by_chr_10/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_R2.rds"))
            #R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/2_multi-lassosum_by_chr_10/Results*_try_truncatedE-3/",ethnic,"_rho_",rho,"_size_",size,"_validation_R2.rds"))
            #R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/2_multi-lassosum_by_chr_10/Results*_try_signif3/",ethnic,"_rho_",rho,"_size_",size,"_validation_R2.rds"))
          }, error=function(e){})
        }else if(method == "5-ethnic lassosum"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/2_multi-lassosum_by_chr_3/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_R2.rds")) }, error=function(e){})
        }else if(method == "2-ethnic lassosum (fix_EUR_coef)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/4_multi-lassosum_by_chr_10_fixEUR/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_R2.rds")) }, error=function(e){})
        }else if(method == "5-ethnic lassosum (fix_EUR_coef)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/4_multi-lassosum_by_chr_3_fixEUR/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_R2.rds")) }, error=function(e){})
        }else if(method == "combined 2-ethnic lassosum"){
          tryCatch({ R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/2_multi-lassosum_by_chr_10/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
        }else if(method == "combined 5-ethnic lassosum"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/2_multi-lassosum_by_chr_3/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
        }else if(method == "combined 2-ethnic lassosum (fix_EUR_coef)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/4_multi-lassosum_by_chr_10_fixEUR/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
        }else if(method == "combined 5-ethnic lassosum (fix_EUR_coef)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/4_multi-lassosum_by_chr_3_fixEUR/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
        }else if(method == "less_para combined 2-ethnic lassosum"){
          tryCatch({ R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/2_multi-lassosum_by_chr_10/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
        }else if(method == "less_para combined 5-ethnic lassosum"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/2_multi-lassosum_by_chr_3/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
        }else if(method == "less_para combined 2-ethnic lassosum (fix_EUR_coef)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/4_multi-lassosum_by_chr_10_fixEUR/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
        }else if(method == "less_para combined 5-ethnic lassosum (fix_EUR_coef)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/4_multi-lassosum_by_chr_3_fixEUR/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
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

res <- data.frame(eth.vec=eth.vec,l_vec=l_vec,m_vec=m_vec,r2.vec=r2.vec,method_vec=method_vec,ga_vec=ga_vec)
res

readr::write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/result_table/Restuls_GA_",setting,"_lassosum2_less_para.txt"))
paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/result_table/Restuls_GA_",setting,"_lassosum2_less_para.txt")



######## new lassosum2

rm(list=ls())

setting=2
rep=1

mulmethods <- c("2-ethnic", "2-ethnic_fixEUR", "5-ethnic", "5-ethnic_fixEUR")
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
      for (method in  c("2-ethnic",
                        "2-ethnic_fixEUR",
                        "5-ethnic",
                        "5-ethnic_fixEUR",

                        "5-ethnic(3)",
                        "5-ethnic(2)",
                        "lassosum2",
                        "EUR lassosum2",
                        "Weighted lassosum2"
                        )){
        i <- i+1
        if(method %in% mulmethods){ mulflag <- T }else{ mulflag <- F }

        if(mulflag){
          if(method == "2-ethnic"){
            tryCatch({R2 <- apply(readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/2_multi-lassosum_by_chr_5/Results/test_all_prs/",ethnic,"/rho_",rho,"_size_",size,"_validation_R2.rds")),MARGIN = 2, mean) }, error=function(e){})
          }else if(method == "2-ethnic_fixEUR"){
            tryCatch({R2 <- apply(readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic,"_EUR/rep_",rep,"/4_multi-lassosum_by_chr_5_fixEUR/Results/test_all_prs/",ethnic,"/rho_",rho,"_size_",size,"_validation_R2.rds")),MARGIN = 2, mean) }, error=function(e){})
          }else if(method == "5-ethnic"){
            tryCatch({R2 <- apply(readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/2_multi-lassosum_by_chr_5/Results/test_all_prs/",ethnic,"/rho_",rho,"_size_",size,"_validation_R2.rds")),MARGIN = 2, mean) }, error=function(e){})
          }else if(method == "5-ethnic_fixEUR"){
            tryCatch({R2 <- apply(readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum_train_from_single_eth/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/4_multi-lassosum_by_chr_5_fixEUR/Results/test_all_prs/",ethnic,"/rho_",rho,"_size_",size,"_validation_R2.rds")),MARGIN = 2, mean) }, error=function(e){})
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
            tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,".rds")) }, error=function(e){})
          }else if(method == "EUR lassosum2"){
            tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_4_EUR.rds")) }, error=function(e){})
          }else if(method == "Weighted lassosum2"){
            tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,"_weighted_lassosum2.rds")) }, error=function(e){})
          }else if(method == "5-ethnic(3)"){
            tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/2_multi-lassosum_by_chr_3/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
          }else if(method == "5-ethnic(2)"){
            tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/AFR_AMR_EAS_EUR_SAS/rep_",rep,"/2_multi-lassosum_by_chr_2/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")) }, error=function(e){})
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
res

readr::write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/result_table/Restuls_GA_",setting,"_lassosum2_from_single_eth.txt"))
paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/result_table/Restuls_GA_",setting,"_lassosum2_from_single_eth.txt")


