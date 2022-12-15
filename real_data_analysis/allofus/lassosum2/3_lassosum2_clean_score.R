# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

NCORES=2

#############################################################
#############################################################

load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/",ethnic,"_",trait,".RData"))

#c(df_beta, beta_lassosum2, params2)

snps <- df_beta$rsid

sumdata <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/",ethnic,"/",trait,".txt"))

sumdataA1 <- sumdata$A1[match(snps, sumdata$rsID)]
sumdataA2 <- sumdata$A2[match(snps, sumdata$rsID)]
m <- df_beta$a1 != sumdataA1

if(sum(m)>0){
  beta_lassosum2[m,] <- -1 * beta_lassosum2[m,]
}

df <- data.frame(SNP = snps, A1= sumdataA1, A2=sumdataA2, beta_lassosum2, stringsAsFactors=F)
fwrite2(df, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/beta_",ethnic,"_",trait,".txt"), col.names = F, sep="\t", nThread=NCORES)
write_tsv(params2, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/para_",ethnic,"_",trait,".txt"))

