# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
M <- length(ethnic)
ethnic1=paste(ethnic,collapse = "_")

### Load snp lists and summary stats
#### Please match summary stats to LD ref alleles before input
Nsnps0 <- vector("list", length = M)
snps_list0 <- vector("list", length = M)
LD_list0 <- vector("list", length = M)
for (l in 1:M){
  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/LDdata/",ethnic[l],"/standard_data/chr",chr,"_LD.RData"))
  Nsnps0[[l]] <- Nsnps
  snps_list0[[l]] <- snps_list
  LD_list0[[l]] <- LD_list
}
nblock <- length(LD_list)
rm(list = c("Nsnps", "snps_list","LD_list"))

## load summary data
summ_list0 <- vector("list", length = M)
summ_max0 <- numeric(length = M)
snps_scale0 <- vector("list", length = M)
N0 <- numeric(length = M)
for (l in 1:M){
  if(ethnic[l]=="EUR"){
    df_beta <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic[l],"/rho_",rho,"_size_",4,"_rep_",rep,"_GA_",setting,"_matched_to_ref.rds"))
  }else{
    df_beta <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic[l],"/rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting,"_matched_to_ref.rds"))
  }
  df_beta$snps_scale <- sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)
  df_beta$beta_hat <- df_beta$beta / df_beta$snps_scale

  N0[l] <- max(df_beta$n_eff)
  summ_max0[l] <- max(abs(df_beta$beta_hat))

  df_beta <- df_beta[df_beta$rsid %in% unlist(snps_list0[[l]]),]
  summ_list0[[l]] <- lapply(snps_list0[[l]], FUN=function (x){ df_beta$beta_hat[match(x, df_beta$rsid)] } )
  snps_scale0[[l]] <- lapply(snps_list0[[l]], FUN=function (x){ df_beta$snps_scale[match(x, df_beta$rsid)] } )

  if(l==1){
    alleles <- df_beta[,c("rsid","a1","a0")]
  }else{
    alleles <- unique(rbind(alleles,df_beta[,c("rsid","a1","a0")]))
  }
  rm(list = c("df_beta"))
}

rm(list = c("df_beta"))

### Load single ethnic best tuning parameter
library(readr)
delta <- numeric(length = M)
lambda <- numeric(length = M)
for (l in 1:M){
  if(ethnic[l]=="EUR"){
    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic[l],"/rep_",rep,"/best_tuning_rho_",rho,"_size_",4,".RData"))
  }else{
    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic[l],"/rep_",rep,"/best_tuning_rho_",rho,"_size_",size,".RData"))
  }
  d <- params2$delta[best]
  tmp <- unique(d)
  tabletmp <- integer(); for (i in 1:length(tmp)) { tabletmp[i] <- sum(d==tmp[i])}
  delta[l] <- tmp[which.max(tabletmp)]

  lam <- params2$lambda[best]
  tmp <- unique(lam)
  tabletmp <- integer(); for (i in 1:length(tmp)) { tabletmp[i] <- sum(lam==tmp[i])}
  lambda[l] <- tmp[which.max(tabletmp)]
}


####################################
## transform into standard data format

indx_block1 <- integer(length = nblock)
snp_list1 <- vector("list", length = nblock)
Nsnps1 <- integer(length = nblock)
indx1 <- vector("list", length = nblock)
summ_list1 <- vector("list", length = nblock)
snps_scale1 <- vector("list", length = nblock)

LD_list1 <- vector("list", length = nblock)
for (bl in 1:nblock){
  ## snp_list1, Nsnps1
  snp_list_tmp <- vector("list", length = M)
  tmp <- character()
  for (l in 1:M){
    if(is.null(snps_list0[[l]][[bl]])) { next }
    snp_list_tmp[[l]] <- snps_list0[[l]][[bl]]
    tmp <- c(tmp, snp_list_tmp[[l]])
  }
  tmp <- unique(tmp)
  Nsnps1[bl] <- length(tmp)

  if(Nsnps1[bl]==0){ indx_block1[bl] <- 0; next }

  snp_list1[[bl]] <- tmp
  indx_block1[bl] <- 1

  ## indx1: the position of ref SNP in original summ_list
  ## summ_list1: summ stat matched to reference snp list (snps not in a certain ethnic group were set to 0)
  ## LD_list1: LD correlations matched to reference snp list (snps not in a certain ethnic group were set to 0)
  indx_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)
  summ_list_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)
  snps_scale_tmp <- matrix(nrow = Nsnps1[bl], ncol = M)

  LD_list_tmp <- vector("list", length = M)
  for (l in 1:M){
    m <- match(snp_list1[[bl]], snp_list_tmp[[l]])
    m1 <- m; m1[is.na(m1)] <- 0
    indx_tmp[,l] <- m1

    m1 <- summ_list0[[l]][[bl]][m]
    summ_list_tmp[,l] <- m1

    m1 <- snps_scale0[[l]][[bl]][m]
    snps_scale_tmp[,l] <- m1
    m1 <- LD_list0[[l]][[bl]][m,m]; m1[is.na(m1)] <- 0; diag(m1)[is.na(m)] <- 1
    LD_list_tmp[[l]] <- m1
  }
  indx1[[bl]] <- indx_tmp
  summ_list1[[bl]] <- summ_list_tmp
  snps_scale1[[bl]] <- snps_scale_tmp
  LD_list1[[bl]] <- LD_list_tmp

  print(bl)
}

summ_list <- summ_list1
snps_scale <- snps_scale1
LD_list <- LD_list1
indx <- indx1
indx_block <- indx_block1
delta_best <- delta
lambda_best <- lambda

snp_list <- snp_list1
Nsnps <- Nsnps1

summ_max <- summ_max0
N <- N0

saveRDS(summ_list,paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/summ_list_rho_",rho,"_size_",size,"_chr",chr,".rds"))
saveRDS(snps_scale,paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/snps_scale_rho_",rho,"_size_",size,"_chr",chr,".rds"))
saveRDS(LD_list,paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/LD_list_rho_",rho,"_size_",size,"_chr",chr,".rds"))
save(indx, indx_block, M, snp_list, alleles, Nsnps, ethnic, delta_best, lambda_best, summ_max, N,
     file=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/otherinfo_rho_",rho,"_size_",size,"_chr",chr,".RData"))

rm(list = c("LD_list","indx","indx_block","Nsnps","delta_best","lambda_best","summ_max","N"))
rm(list = c("LD_list1","indx1","indx_block1","Nsnps1","summ_max0","N0"))

################################################################
################################################################

## fixEUR data

summ_list <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/summ_list_rho_",rho,"_size_",size,"_chr",chr,".rds"))
snps_scale <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/snps_scale_rho_",rho,"_size_",size,"_chr",chr,".rds"))
load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/otherinfo_rho_",rho,"_size_",size,"_chr",chr,".RData"))

chr <- as.integer(chr)

### Load EUR coef in single ethnic lassosum2
# df_beta, beta_lassosum2, params2,
load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/EUR/rep_",rep,"/beta_rho_",rho,"_size_4.RData"))
m <- df_beta$chr==chr
df_beta <- df_beta[m,]; beta_lassosum2 <- beta_lassosum2[m,]

library(dplyr)
ref <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/EUR/ref_chr",chr,".bim"))
snpinfo <- left_join(df_beta[,c(5,4)], ref[,c(2,5)], by=c("rsid"="V2"))
m <- snpinfo$a1!=snpinfo$V5; if(sum(m)>0){ beta_lassosum2[m,] <- -beta_lassosum2[m,] }

### Load best tuning
load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/EUR/rep_",rep,"/best_tuning_rho_",rho,"_size_4.RData"))
tmp <- unique(best); tabletmp <- integer(); for (i in 1:length(tmp)) { tabletmp[i] <- sum(best==tmp[i])}; best <- tmp[which.max(tabletmp)]
b_EUR <- beta_lassosum2[,best]

### convert original scale effect sizes output from lassosum2 to standardized scale
m <- which(ethnic=="EUR")
for (i in 1:length(snp_list)){
  if(length(snp_list[[i]])>0){
    tmp <- b_EUR[match(snp_list[[i]],df_beta$rsid)]
    summ_list[[i]][,m] <- tmp/snps_scale[[i]][,m]
  }
  print(i)
}

saveRDS(summ_list,paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/summ_list_fixEUR_rho_",rho,"_size_",size,"_chr",chr,".rds"))

