
rm(list=ls())

setting=5; rep=1
ethnic_folder='AFR_AMR_EAS_EUR_SAS'
r2.vec <- numeric(); eth.vec <- character(); l_vec <- integer(); m_vec <- integer(); method_vec <- character(); ga_vec <- integer()
i <- 0
for (ethnic in c("AFR","AMR","EAS","SAS")){
  for (method in  c("multi-ethnic lasso (freelambda_10_10)",
                    "multi-ethnic lasso (fix_EUR_beta)",
                    "multi-ethnic lasso (all-ethnic_3)",
                    "multi-ethnic lasso^2 (freelambda_10_10)",
                    "multi-ethnic lasso^2 (fix_EUR_beta)",
                    "multi-ethnic lasso^2 (all-ethnic_3)",
                    "lassosum",
                    "lassosum2",
                    "EUR lassosum",
                    "EUR lassosum2",
                    "Weighted-PRS (lassosum)",
                    "Weighted-PRS (lassosum2)")){
    for (rho in 1:3){
      for (size in 1:4){
        i <- i+1
        if(method == "multi-ethnic lasso (freelambda_10_10)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"/RunPRS/2_MultiEthnic/", ethnic, "/freelambda_2ethnicity_10_10_rep_",rep,"/Results/rho_",rho,"_size_",size,"_validation_R2.rds")) }, error=function(e){})
        }else if(method == "multi-ethnic lasso (fix_EUR_beta)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"/RunPRS/2_MultiEthnic/", ethnic, "/fixEUR_2ethnicity_rep_",rep,"/Results/rho_",rho,"_size_",size,"_validation_R2.rds")) }, error=function(e){})
        }else if(method == "multi-ethnic lasso (all-ethnic_3)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic/Setting_",setting,"/RunPRS/2_MultiEthnic/",ethnic_folder,"/rep_",rep,"/5ethnicity_3/Results/",ethnic,"_rho_",rho,"_size_",size,"_validation_R2.rds")) }, error=function(e){})
        }else if(method == "multi-ethnic lasso^2 (freelambda_10_10)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"/RunPRS/2_MultiEthnic/",ethnic, "/freelambda_2ethnicity_10_10_rep_",rep,"/Results/rho_",rho,"_size_",size,"_validation_lassoSCORE_R2.rds")) }, error=function(e){})
        }else if(method == "multi-ethnic lasso^2 (fix_EUR_beta)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"/RunPRS/2_MultiEthnic/",ethnic, "/fixEUR_2ethnicity_rep_",rep,"/Results/rho_",rho,"_size_",size,"_validation_lassoSCORE_R2.rds")) }, error=function(e){})
        }else if(method == "multi-ethnic lasso^2 (all-ethnic_3)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic/Setting_",setting,"/RunPRS/2_MultiEthnic/",ethnic_folder,"/rep_",rep,"/5ethnicity_3/Results/",ethnic2,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2.rds")) }, error=function(e){})
        }else if(method == "lassosum"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"/RunPRS/1_SingleEthnic/",ethnic, "/", ethnic,"/rep_",rep,"/Results/rho_",rho,"_size_",size,"_validation_R2.rds")) }, error=function(e){})
        }else if(method == "lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,".rds")) }, error=function(e){})
        }else if(method == "EUR lassosum"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"/RunPRS/1_SingleEthnic/",ethnic, "/", ethnic,"/rep_",rep,"/Results/rho_",rho,"_size_",size,"_validation_EUR_R2.rds")) }, error=function(e){})
        }else if(method == "EUR lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_4_EUR.rds")) }, error=function(e){})
        }else if(method == "Weighted-PRS (lassosum)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"/RunPRS/1_SingleEthnic/",ethnic, "/",ethnic, "/rep_",rep,"/Results/rho_",rho,"_size_",size,"_validation_weighted_lassoPRS_R2.rds")) }, error=function(e){})
        }else if(method == "Weighted-PRS (lassosum2)"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,"_weighted_lassosum2.rds")) }, error=function(e){})
        }

        tryCatch({r2.vec[i] <- mean(R2) }, error=function(e){})
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
readr::write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/result_table/Restuls_GA_",setting,"_lassosum2.txt"))
paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/result_table/Restuls_GA_",setting,"_lassosum2.txt")
