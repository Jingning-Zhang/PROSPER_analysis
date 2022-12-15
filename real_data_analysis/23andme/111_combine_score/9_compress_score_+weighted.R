# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)
library(dplyr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

#NCORES <- 3
#digits <- 3


if(targetethnic=="EUR"){
  targetethnic1 <- "european"
}else if(targetethnic=="AFR"){
  targetethnic1 <- "african_american"
}else if(targetethnic=="AMR"){
  targetethnic1 <- "latino"
}else if(targetethnic=="EAS"){
  targetethnic1 <- "east_asian"
}else if(targetethnic=="SAS"){
  targetethnic1 <- "south_asian"
}

#############################################################
#############################################################

## binded_file: 5eth m2 (multi-lassosum_train_from_single_eth) + weighted + single test + baseline test

dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/"))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/binded_file/"))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/binded_file/",targetethnic))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/binded_file/",targetethnic,"/",trait))

dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/compressed_binded_file/"))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/compressed_binded_file/"))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/compressed_binded_file/",targetethnic1))

########################
## 5eth m2
#methods = c("multi-lassosum/Results_to_23andme/m1_2ethnicity",
#            "multi-lassosum_train_from_single_eth/Results_to_23andme/m2_2ethnicity",
#            "multi-lassosum_train_from_single_eth/Results_to_23andme/m2_5ethnicity")
#methods1 = c("m1_2ethnicity",
#            "m2_2ethnicity",
#            "m2_5ethnicity")
#for (i in 1:3){
#  print(i)
  df <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/multi-lassosum_train_from_single_eth/Results_to_23andme/m2_5ethnicity/",
                                targetethnic,"/",trait,"/prs.file"))
  colnames(df)[-1] <- paste0("m2_5ethnicity-",1:(ncol(df)-1))
  #if(i==1){
      dffull <- df
  #}else{
  #  dffull <- full_join(dffull, df, by=c("V1"))
  #}
#}

########################
## weighted

for (weight_ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  df <- fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_best/",weight_ethnic,"_",trait,".prs.file"),
               sep=" ")
  if(weight_ethnic == "EUR"){
    weight_dffull <- df
  }else{
    weight_dffull <- full_join(weight_dffull, df, by=c("V1"))
  }
}
colnames(weight_dffull)[-1] <- paste0("weighted-",1:(ncol(weight_dffull)-1))

dffull <- full_join(dffull, weight_dffull, by=c("V1"))

########################
## single test

df <- fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_best/",targetethnic,"_",trait,".prs.file"),
               sep=" ")

singletest_dffull <- df

colnames(singletest_dffull)[-1] <- paste0("single-1")

dffull <- full_join(dffull, singletest_dffull, by=c("V1"))

########################
## baseline test

df <- fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_best/","EUR","_",trait,".prs.file"),
               sep=" ")

baseline_dffull <- df

colnames(baseline_dffull)[-1] <- paste0("baseline-1")

dffull <- full_join(dffull, baseline_dffull, by=c("V1"))


########################
########################
## final clean

for (i in 2:ncol(dffull)){
  dffull[[i]][is.na(dffull[[i]])] <- 0
  dffull[[i]][(dffull[[i]] == Inf) | (dffull[[i]] == -Inf)] <- 0
}
dffull[[1]] <- format(dffull[[1]], trim = T, scientific = F)

df_order <- colnames(dffull)
fwrite2(dffull, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/binded_file/",targetethnic,"/",trait,"/prs.file"),
        col.names = F, sep=" ", nThread=1)
writeLines(df_order, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/binded_file/",targetethnic,"/",trait,"/prs.order"))


setwd(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/binded_file/",targetethnic))
system(paste0("tar -czvf ",trait,".tar.gz ",trait))
system(paste0("mv ",trait,".tar.gz", " /dcs04/nilanjan/data/23andme/Analysis/JZ/Results_to_23andme_2022_Nov/compressed_binded_file/",targetethnic1,"/"))



