

rm(list=ls())

library(readr)

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/clean_results")

r2.vec <- numeric(); eth.vec <- character(); t_vec <- character(); method_vec <- character()
i <- 0

#m1="multi-lassosum"
m1="multi-lassosum_train_from_single_eth"
#m1="multi-lassosum_two_parameter"

if(m1=="multi-lassosum"){ para1=3 }else{ para1=10 }

path = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus",
              "/",m1,
              "/AFR_AMR_EUR",
              "/2_multi-lassosum_by_chr_",para1,
              "/Results/final_results/")

for (ethnic in c("AFR","AMR")){
  for (method in  c("lassosum2",
                    "EUR lassosum2",
                    "weighted lassosum2",
                    "multi-ethnic lasso")){
    for (trait in c("height","bmi")){

        if(method == "lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/final_results/",ethnic,"_",trait,"_validation_R2_best.rds")) }, error=function(e){})

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,1] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- method

        }else if(method == "EUR lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/final_results/",ethnic,"_",trait,"_validation_R2_best_EUR_lassosum2.rds")) }, error=function(e){})

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,1] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- method

        }else if(method == "weighted lassosum2"){
          tryCatch({R2 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/weighted_lassosum2/final_results/",ethnic,"_",trait,"_validation_R2_weighted.rds")) }, error=function(e){})

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,1] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (2)")

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,2] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (all ancestry)")
        }else if(method == "multi-ethnic lasso"){
          tryCatch({R2 <- readRDS(paste0(path,"/",ethnic,"/",trait,"_validation_R2.rds")) }, error=function(e){})

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,2] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (best)")

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,4] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (lasso)")

          i <- i+1
          a <- tryCatch({r2.vec[i] <- R2[1,6] }, error=function(e){})
          if(is.null(a)){r2.vec[i] <- NA}
          eth.vec[i] <- ethnic
          t_vec[i] <- trait
          method_vec[i] <- paste0(method, " (super learning)")
        }

        rm(R2)

    }
  }
}

res <- data.frame(eth.vec=eth.vec,t_vec=t_vec,r2.vec=r2.vec,method_vec=method_vec)

readr::write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/clean_results/table_",m1,".txt"))
paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/clean_results/table_",m1,".txt")

