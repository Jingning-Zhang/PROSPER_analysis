# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

NCORES <- 1


library(dplyr)
for (trait in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){
  for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
    df <- fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_best/",ethnic,"_",trait,".prs.file"),
                 sep=" ",
                 colClasses = c("character",  "numeric" ))
    colnames(df) <- c("V1",ethnic)
    if(ethnic=="EUR"){
      dffull <- df
    }else{
      dffull <- full_join(dffull, df, by="V1")
    }
  }
  dffull$EUR[is.na(dffull$EUR)] <- 0
  dffull$AFR[is.na(dffull$AFR)] <- 0
  dffull$AMR[is.na(dffull$AMR)] <- 0
  dffull$EAS[is.na(dffull$EAS)] <- 0
  dffull$SAS[is.na(dffull$SAS)] <- 0
  fwrite2(dffull, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/weighted_lassosum2/all_ethnic_",trait,".weighted.prs.file"), col.names = F, sep=" ", nThread=1)

}


#############################################################
## match to 23andme SNP ID

load("/dcs04/nilanjan/data/23andme/cleaned/snpinfo/snpinfo_mega.RData")

ethnic='AFR'
trait='height'

sumdata <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/23andme/cleaned/",ethnic,"/sumdat/",trait,"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"))
sumdata <- sumdata[,c(1,4)]

############
# two ethnic
wprs <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/weighted_lassosum2/",ethnic,"_",trait,".weighted.prs.file"))
rsid <- snpinfo_mega$assay.name[match(wprs$V1,snpinfo_mega$im.data.id)]
A1 <- sumdata$A1[match(rsid,sumdata$rsid)]
res <- data.frame(rsid, A1, wprs[,2:3])
fwrite2(res, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/weighted_lassosum2/ukbb_test/",ethnic,"_",trait,".prs"), col.names = F, sep=" ", nThread=NCORES)


############
# all ethnic
wprs <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/weighted_lassosum2/all_ethnic_",trait,".weighted.prs.file"))
rsid <- snpinfo_mega$assay.name[match(wprs$V1,snpinfo_mega$im.data.id)]
A1 <- sumdata$A1[match(rsid,sumdata$rsid)]
res <- data.frame(rsid, A1, wprs[,2:6])
fwrite2(res, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/weighted_lassosum2/ukbb_test/all_ethnic_",trait,".prs"), col.names = F, sep=" ", nThread=NCORES)



