# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

NCORES <- 1

#############################################################
#############################################################

M <- length(ethnic)
ethnic1 <- paste(ethnic,collapse = "_")

whichethnic <- which(ethnic == targetethnic)


df_beta <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/",targetethnic,"_",trait,"_df_beta.rds"))
snps_scale <- data.frame(rsid=df_beta$rsid, snps_scale=sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2))

for (chr in 1:22){
  print(chr)
  # ethnic, M, delta_best, runtime, res,
  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/tuning/",trait,"_chr",chr,".RData"))
  delta <- delta_best
  #indx, indx_block, M, snp_list, ethnic,
  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/chr",chr,".RData"))

  if(chr==1){
    alltuning <- length(res$c)
    conv <- matrix(nrow=alltuning,ncol = 22)

    ## param
    param <- matrix(nrow=alltuning, ncol = 2*M+2)
    for (m in 1:M){ param[,m] <- delta[m] }
    param[,(M+1):(2*M)] <- t(res$lambda)
    param[,(2*M+1)] <- unlist(lapply(res$c, FUN = function (x){x[1,2]}))

    ## df_beta_allchr
    snps <- unlist(snp_list)
    m_snps <- snps %in% snps_scale$rsid
    bim.ref <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/ref_geno/",targetethnic,"/",trait,"/ref_chr", chr,".bim"))
    bim.ref <- bim.ref[match(snps, bim.ref$V2),]
    A1 <- bim.ref$V5
    df_beta_allchr <- data.frame(rsid=snps[m_snps], a1=A1[m_snps], stringsAsFactors = F)

    ## b_allchr
    b_tmp <- matrix(nrow = sum(m_snps), ncol = alltuning)
    for (i in 1:alltuning){ b_tmp[,i] <- unlist(lapply(res$b[[i]], FUN=function (x) {x[,whichethnic]}))[m_snps] * snps_scale$snps_scale[match(snps[m_snps],snps_scale$rsid)] } # dosage effect
    b_allchr <- b_tmp

    ## conv
    res$conv[res$conv==-2147483648] <- 1
    conv[,chr] <- as.integer(apply(res$conv,MARGIN = 2, FUN = function (x){mean(x==1)==1}))

  }else{

    ## df_beta_allchr
    snps <- unlist(snp_list)
    m_snps <- snps %in% snps_scale$rsid
    bim.ref <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/ref_geno/",targetethnic,"/",trait,"/ref_chr", chr,".bim"))
    bim.ref <- bim.ref[match(snps, bim.ref$V2),]
    A1 <- bim.ref$V5
    df_beta_allchr <- rbind(df_beta_allchr,
                            data.frame(rsid=snps[m_snps], a1=A1[m_snps], stringsAsFactors = F))

    ## b_allchr
    b_tmp <- matrix(nrow = sum(m_snps), ncol = alltuning)
    for (i in 1:alltuning){ b_tmp[,i] <- unlist(lapply(res$b[[i]], FUN=function (x) {x[,whichethnic]}))[m_snps] * snps_scale$snps_scale[match(snps[m_snps],snps_scale$rsid)] } # dosage effect
    b_allchr <- rbind(b_allchr, b_tmp)

    ## conv
    res$conv[res$conv==-2147483648] <- 1
    conv[,chr] <- as.integer(apply(res$conv,MARGIN = 2, FUN = function (x){mean(x==1)==1}))
  }

  rm(list = c("bim.ref", "b_tmp", "res","indx","indx_block","snp_list"))

}
param[,(2*M+2)] <- apply(b_allchr, MARGIN = 2, FUN = function (x){mean(x!=0)})
colnames(param) <- c(paste0("delta",1:M), paste0("lambda",1:M), "c", "sparsity")

conv_allchr <- apply(conv, MARGIN = 1, sum); m <- conv_allchr == 22
if(sum(!m)!=0){
  param <- param[m,]
  b_allchr <- b_allchr[,m]
}


#############################################
## format results

snps <- df_beta_allchr$rsid

sumdata <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/",ethnic,"/",trait,".txt"))
sumdataA1 <- sumdata$A1[match(snps, sumdata$rsID)]
sumdataA2 <- sumdata$A2[match(snps, sumdata$rsID)]
m <- df_beta_allchr$a1 != sumdataA1
if(sum(m)>0){ b_allchr[m,] <- -1 * b_allchr[m,] }

df <- data.frame(SNP = snps, A1= sumdataA1, A2= sumdataA2, b_allchr, stringsAsFactors=F)
fwrite2(df, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/Results/beta_",targetethnic,"_",trait,".txt"), col.names = F, sep="\t", nThread=NCORES)
write_tsv(as.data.frame(param), paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/Results/para_",trait,".txt"))

rm(list = c("df","param","sumdata"))

