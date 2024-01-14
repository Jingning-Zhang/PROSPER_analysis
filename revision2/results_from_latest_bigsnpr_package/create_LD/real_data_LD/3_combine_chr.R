
rm(list=ls())
library(bigsnpr)
library(bigreadr)

races = c("AFR","AMR","EAS","EUR","SAS")

#args <- commandArgs(T)
#for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
#
#race = races[as.numeric(temp1)]


for(race in races){

  path = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/",race)
  setwd(path)

  map_ldref <- NULL
  for(chr in 1:22){
    print(paste0(race,"-",chr))
    map_ldref <- rbind(map_ldref, readRDS(paste0("./map_ldref_chr",chr,".rds")))
  }
  saveRDS(map_ldref, paste0("./map_ldref.rds"))
}

