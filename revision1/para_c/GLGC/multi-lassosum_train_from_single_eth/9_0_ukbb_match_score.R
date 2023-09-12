# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

#ethnic="AFR"
#trait="HDL"

rm(list=ls())

allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)

  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;

  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  snp[["amb_str"]] = (a1 == flip2 & a2 == flip1) | (a1 == flip1 & a2 == flip2)

  return(snp)
}

library(readr)
library(bigreadr)


args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
L <- as.numeric(L); Lc <- as.numeric(Lc)

NCORES <- 1

#############################################################
#############################################################

ethnic1 <- paste(ethnic,collapse = "_")

path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/para_c/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para)

#############################################################
#############################################################

df <- fread2(paste0(path,"/Results/beta_",targetethnic,"_",trait,".txt"), header = F, sep="\t", nThread=NCORES)
snps <- df$V1
sumdataA1 <- df$V2
sumdataA2 <- df$V3
beta_lassosum2 <- as.matrix(df[,-1:-3])

bim <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/",ethnic,"/allchr_bim"))

m <- which(!(snps %in% bim$V2))
bim <- bim[match(snps, bim$V2),]

if(length(m)>0) {
  snps <- snps[-m]
  sumdataA1 <- sumdataA1[-m]
  sumdataA2 <- sumdataA2[-m]
  bim <- bim[-m,]
  beta_lassosum2 <- beta_lassosum2[-m,]
}

qc <- allele.qc(a1=sumdataA1,a2=sumdataA2,ref1=bim$V5,ref2=bim$V6)

if(sum(qc$flip)>0){
  beta_lassosum2[qc$flip,] <- -1 * beta_lassosum2[qc$flip,]
}

df <- data.frame(SNP = snps, A1= bim$V5, beta_lassosum2, stringsAsFactors=F)
fwrite2(df, paste0(path,"/Results/ukbb_beta_",targetethnic,"_",trait,".txt"), col.names = F, sep="\t", nThread=NCORES)

alltuning <- ncol(beta_lassosum2)
saveRDS(alltuning, paste0(path,"/Results/alltuning_",targetethnic,"_",trait,".rds"))
