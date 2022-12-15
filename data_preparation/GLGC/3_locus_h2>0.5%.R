

windowsize=500000

library(dplyr)
library(bigreadr)
library(readr)
trait='HDL'
ethnic='EUR'
for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
  df_hugeh2 <- tibble()
  for(ethnic in c("EUR","AFR","AMR","EAS","SAS")){
    print(ethnic)
    sumdata <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",ethnic,"/",trait,".txt"))
    df_beta = sumdata
    df_beta$snps_scale <- sqrt(df_beta$N * df_beta$SE^2 + df_beta$BETA^2)
    df_beta$beta_hat <- df_beta$BETA / df_beta$snps_scale
    df_beta$Z <- df_beta$BETA / df_beta$SE
    tmp <- df_beta[abs(df_beta$beta_hat)^2>0.005,]
    tmp$ancestry <- ethnic
    df_hugeh2 <- rbind(df_hugeh2, tmp)
  }
  df_hugeh2 <- df_hugeh2[,c("rsID","CHR","POS_b37","N","P","beta_hat","Z","ancestry")]
  write_tsv(df_hugeh2, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,".txt"))
}


for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
  df_hugeh2 <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,".txt"), col_type=cols())
  loc_chr <- unique(df_hugeh2$CHR)
  res <- data.frame(chr=0,start=0,end=0,topSNP="a",z=0,ancestry="a")
  for(i in 1:length(loc_chr)){
    df_hugeh2_tmp <- df_hugeh2[df_hugeh2$CHR==loc_chr[i],]
    flag=T
    while(flag){
      m <- which.max(abs(df_hugeh2_tmp$Z))
      res <- rbind(res,
                 data.frame(chr=loc_chr[i],
                            start=df_hugeh2_tmp$POS_b37[m] - windowsize,
                            end=df_hugeh2_tmp$POS_b37[m] + windowsize,
                            topSNP=df_hugeh2_tmp$rsID[m],
                            z=df_hugeh2_tmp$Z[m],
                            ancestry=df_hugeh2_tmp$ancestry[m]))
      mm <- (df_hugeh2_tmp$POS_b37 < df_hugeh2_tmp$POS_b37[m] - windowsize) & (df_hugeh2_tmp$POS_b37 > df_hugeh2_tmp$POS_b37[m] + windowsize)
      if(sum(mm)==0){
        flag=F
      }
    }
  }
  res <- res[-1,]
  write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,"_topSNP.txt"))
}

for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
  res <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,"_topSNP.txt"), col_type=cols())
  sumdata <- list()
  for(ethnic in c("EUR","AFR","AMR","EAS","SAS")){
      sumdata[[ethnic]] <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",ethnic,"/",trait,".txt"))
  }
  for(i in 1:nrow(res)){
    flag=T
    for(ethnic in c("EUR","AFR","AMR","EAS","SAS")){
      if(!(res$topSNP[i] %in% sumdata[[ethnic]]$rsID)){flag=F;break}
    }
    if(!flag){
      print(paste0(trait,":",i))
      snps <- list()
      for(ethnic in c("EUR","AFR","AMR","EAS","SAS")){
        snps[[ethnic]] <- sumdata[[ethnic]]$rsID[(sumdata[[ethnic]]$CHR==res$chr[i]) & (sumdata[[ethnic]]$POS_b37 >= res$start[i]) & (sumdata[[ethnic]]$POS_b37 <= res$end[i])]
      }
      intersnps <- snps[[1]]
      for(j in 2:length(snps)){
        intersnps <- intersect(intersnps, snps[[j]])
      }
      sumdata_ref <- sumdata[[res$ancestry[i]]][sumdata[[res$ancestry[i]]]$rsID %in% intersnps,]
      sumdata_ref$Z <- sumdata_ref$BETA / sumdata_ref$SE
      m <- which.max(abs(sumdata_ref$Z))
      res$start[i] <- sumdata_ref$POS_b37[m] - windowsize
      res$end[i] <- sumdata_ref$POS_b37[m] + windowsize
      res$topSNP[i] <- sumdata_ref$rsID[m]
      res$z[i] <- sumdata_ref$Z[m]
    }
  }
  colnames(res)[4] <- "snp"
  write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,"_topSNP_corrected.txt"))
}

for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
  res <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/",trait,"_topSNP_corrected.txt"), col_type=cols())

  for(ethnic in c("EUR","AFR","AMR","EAS","SAS")){
    sumdata <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",ethnic,"/",trait,".txt"))

    for(i in 1:nrow(res)){
      if(i==1){
        mm <- (sumdata$CHR==res$chr[i]) & (sumdata$POS_b37 >= res$start[i]) & (sumdata$POS_b37 <= res$end[i])
      }else{
        mm <- mm | ( (sumdata$CHR==res$chr[i]) & (sumdata$POS_b37 >= res$start[i]) & (sumdata$POS_b37 <= res$end[i]) )
      }
    }
    sumdata_remove <- sumdata[!mm,]

    write_tsv(sumdata_remove, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",ethnic,"/",trait,"_remove_big_effect_locus.txt"))

  }

}


####### biggest snp

for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
  df_hugeh2 <- tibble()
  for(ethnic in c("EUR","AFR","AMR","EAS","SAS")){
    print(ethnic)
    sumdata <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",ethnic,"/",trait,".txt"))
    df_beta = sumdata
    df_beta$snps_scale <- sqrt(df_beta$N * df_beta$SE^2 + df_beta$BETA^2)
    df_beta$beta_hat <- df_beta$BETA / df_beta$snps_scale
    df_beta$Z <- df_beta$BETA / df_beta$SE
    tmp <- df_beta[which.max((df_beta$beta_hat)^2),]
    tmp$ancestry <- ethnic
    tmp$h2 <- (tmp$beta_hat)^2
    df_hugeh2 <- rbind(df_hugeh2, tmp)
  }
  df_hugeh2 <- df_hugeh2[,c("rsID","CHR","POS_b37","N","P","beta_hat","Z","ancestry","h2")]
  df_hugeh2 <- df_hugeh2[,c("rsID","ancestry","h2")]
  write_tsv(df_hugeh2, paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/hugeh2_locus/top_h2_",trait,".txt"))
}
