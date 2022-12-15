# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

NCORES <- 1

#trait='height'
#ethnic=c("AFR","EUR")
#L='4'
#Lc='4'
#para='4'
#targetethnic='AFR'


M <- length(ethnic)
ethnic1 <- paste(ethnic,collapse = "_")


load("/dcs04/nilanjan/data/23andme/cleaned/snpinfo/snpinfo_mega.RData")

library(dplyr)

  for (i in 1:length(ethnic)){
    df <- fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/Results/beta_",trait,"_",ethnic[i],".txt"),
                       sep="\t")
    colnames(df) <- c("ID","A1",paste0(ethnic[i],"_",1:(ncol(df)-2)) )

    if(i==1){
      dffull <- df
    }else{
      dffull <- full_join(dffull, df)
    }
  }

  for (i in 3:ncol(dffull)){
    dffull[[i]][is.na(dffull[[i]])] <- 0
  }

ID <- dffull[[1]]
A1 <- dffull[[2]]
dffull <- as.matrix(dffull[,-1:-2])

tmp <- apply(dffull, MARGIN=1, function(x){sum(x!=0)}); m <- !(tmp==0)

ID <- ID[m]
A1 <- A1[m]
dffull <- dffull[m,,drop=F]

id_23andme <- format(snpinfo_mega$im.data.id[match(ID, snpinfo_mega$assay.name)], trim = T, scientific = F)
  rm(list = c("snpinfo_mega","df"))

df <- data.frame(SNP = id_23andme, signif(dffull, digits = 3))
fwrite2(df, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/multi-lassosum_train_from_single_eth/Results_to_23andme/m2_",M,"ethnicity/",targetethnic,"/",trait,"/prs.file"),
        col.names = F, sep=" ", nThread=1)



