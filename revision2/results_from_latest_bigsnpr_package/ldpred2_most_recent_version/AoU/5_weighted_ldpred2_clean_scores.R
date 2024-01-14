# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

NCORES <- 1

library(dplyr)
for (trait in c('height', 'bmi')){
  setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait))

  for (race in c("EUR","AFR","AMR")){
    print(paste0(trait,"-",race))
    for(chr in 1:22){
      tmp <- fread2(paste0("./ldpred2/coef/",race,"_ldpred2-chr",chr,".txt"))
      best <- readRDS(paste0("./ldpred2/prs/tuning/",race,"_best.rds"))
      tmp <- tmp[,c(1,2,best+2)]
      colnames(tmp) <- c("rsID","A1",race)
      if(chr==1){ df <- tmp }else{ df <- rbind(df, tmp) }
    }
    if(race=="EUR"){
      dffull <- df
    }else{
      dffull <- full_join(dffull, df, by=c("rsID","A1"))
    }
    rm(list=c("df","best"))
  }
  for(i in 3:ncol(dffull)){
    dffull[[i]][is.na(dffull[[i]])] <- 0
  }
  fwrite2(dffull, paste0("./ldpred2/coef/all_ancs_best_ldpred2.txt"), sep=" ", nThread=1)
  rm(list="dffull")
}

