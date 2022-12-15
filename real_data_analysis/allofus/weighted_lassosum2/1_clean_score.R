# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

NCORES <- 1

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/weighted_lassosum2")
library(dplyr)
for (trait in c("height","bmi")){
  for (ethnic in c("EUR","AFR","AMR")){
    print(paste0(trait,"-",ethnic))
    df <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/ukbb_beta_",ethnic,"_",trait,".txt"))
    best <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/final_results/",ethnic,"_",trait,"_best.rds"))
    df <- df[,c(1,2,best+2)]
    colnames(df) <- c("rsID","A1",ethnic)
    if(ethnic=="EUR"){
      dffull <- df
    }else{
      dffull <- full_join(dffull, df, by=c("rsID","A1"))
    }
    rm(list=c("df","best"))
  }
  for(i in 3:ncol(dffull)){
    dffull[[i]][is.na(dffull[[i]])] <- 0
  }
  fwrite2(dffull, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/weighted_lassosum2/ukbb_beta_",trait,".txt"), sep=" ", nThread=1)
  rm(list="dffull")
}

