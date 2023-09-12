

rm(list=ls())

library(readr)

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/clean_results")

r2.vec <- numeric(); eth.vec <- character(); t_vec <- character(); method_vec <- character()
i <- 0

i=0
for (ethnic in c("AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","TC")){
    i <- i+1
    a <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC",
              "/multi-lassosum_train_from_single_eth",
              "/AFR_AMR_EAS_SAS_EUR",
              "/2_multi-lassosum_by_chr_10",
              "/Results/final_results/",ethnic,"/",trait,"_validation_R2.rds"))[["superlearning_combined_R2_cross_ancestry"]]
    r2.vec[i] <- a
    eth.vec[i] <- ethnic
    t_vec[i] <- trait
    method_vec[i] <- "original"

    i <- i+1
    b <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/para_c/GLGC",
            "/multi-lassosum_train_from_single_eth/",
            "/AFR_AMR_EAS_SAS_EUR",
            "/2_multi-lassosum_by_chr_10",
            "/Results/final_results/",ethnic,"/",trait,"_validation_R2.rds"))[["superlearning_combined_R2_cross_ancestry"]]
    r2.vec[i] <- b
    eth.vec[i] <- ethnic
    t_vec[i] <- trait
    method_vec[i] <- "flexible_c"

  }
}

res <- data.frame(eth.vec=eth.vec,t_vec=t_vec,r2.vec=r2.vec,method_vec=method_vec)

readr::write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/clean_results/table_",m1,".txt"))
paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/clean_results/table_",m1,".txt")

