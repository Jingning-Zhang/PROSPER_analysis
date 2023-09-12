

rm(list=ls())

library(readr)


r2.vec <- numeric(); eth.vec <- character(); t_vec <- character(); method_vec <- character()
i <- 0

#m1="multi-lassosum"
m1="multi-lassosum_train_from_single_eth"
#m1="multi-lassosum_two_parameter"

if(m1=="multi-lassosum"){ para1=2 }else{ para1=10 }

path = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC",
              "/",m1,
              "/AFR_AMR_EAS_SAS_EUR",
              "/2_multi-lassosum_by_chr_",para1,
              "/Results/final_results/")

for (ethnic in c("AFR","AMR","EAS","SAS")){
  for (method in  c("lassosum2",
                    "PROSPER")){
    for (trait in c("HDL","LDL","logTG","nonHDL","TC")){

        if(method == "lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/final_results/",ethnic,"_",trait,"_validation_R2_sl.rds")) }, error=function(e){})

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,1] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (best)")

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,3] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (sl)")

        }else if(method == "PROSPER"){
          tryCatch({R2 <- readRDS(paste0(path,"/",ethnic,"/",trait,"_validation_R2.rds")) }, error=function(e){})

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,2] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (best)")

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,6] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (sl)")
        }

        rm(R2)

    }
  }
}

res <- data.frame(eth.vec=eth.vec,t_vec=t_vec,r2.vec=r2.vec,method_vec=method_vec)

readr::write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/lassosum2_sl/GLGC_lassosum2_sl.txt"))
paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/lassosum2_sl/GLGC_lassosum2_sl.txt")

