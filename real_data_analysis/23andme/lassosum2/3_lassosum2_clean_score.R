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

load(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/",ethnic,"_",trait,".RData"))

#c(df_beta, beta_lassosum2, params2)

snps <- df_beta$rsid

sumdata <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/23andme/cleaned/",ethnic,"/sumdat/",trait,"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"))
sumdata <- sumdata[,c(1,4)]

sumdataA1 <- sumdata$A1[match(snps, sumdata$rsid)]
m <- df_beta$a1 != sumdataA1

if(sum(m)>0){
  beta_lassosum2[m,] <- -1 * beta_lassosum2[m,]
}
df <- data.frame(SNP = snps, A1= sumdataA1, beta_lassosum2)
fwrite2(df, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/beta_",ethnic,"_",trait,".txt"), col.names = F, sep="\t", nThread=NCORES)
write_tsv(params2, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/para_",ethnic,"_",trait,".txt"))

rm(list = c("df", "df_beta","params2","sumdata"))

#############################################################
## match to 23andme SNP ID

load("/dcs04/nilanjan/data/23andme/cleaned/snpinfo/snpinfo_mega.RData")
id_23andme <- format(snpinfo_mega$im.data.id[match(snps,snpinfo_mega$assay.name)], trim = T, scientific = F)

tmp <- apply(beta_lassosum2, MARGIN=1, function(x){sum(x!=0)})
m <- !(tmp==0)
beta_lassosum2 <- beta_lassosum2[m,]
id_23andme <- id_23andme[m]
rm(list = c("snpinfo_mega"))
df <- data.frame(SNP = id_23andme, beta_lassosum2)
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results_to_23andme/"))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results_to_23andme/",ethnic))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results_to_23andme/",ethnic,"/",trait))
fwrite2(df, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results_to_23andme/",ethnic,"/",trait,"/prs.file"), col.names = F, sep=" ", nThread=NCORES)

