# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
M <- length(ethnic)
ethnic1=paste(ethnic,collapse = "_")

### Load snp lists and summary stats
Nsnps0 <- vector("list", length = M)
snps_list0 <- vector("list", length = M)
LD_list0 <- vector("list", length = M)
for (l in 1:M){
  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic[l],"/",trait,"/standard_data/chr",chr,"_LD.RData"))
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
  df_beta <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/",ethnic[l],"_",trait,"_df_beta.rds"))
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
  para <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/para_",ethnic[l],"_",trait,".txt"))
  best <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/final_results/",ethnic[l],"_",trait,"_best.rds"))
  delta[l] <- para$delta[best]
  lambda[l] <- para$lambda[best]
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
    m1 <- summ_list0[[l]][[bl]][m]; m1[is.na(m1)] <- 0
    summ_list_tmp[,l] <- m1
    m1 <- snps_scale0[[l]][[bl]][m]; m1[is.na(m1)] <- 0
    snps_scale_tmp[,l] <- m1
    m1 <- LD_list0[[l]][[bl]][m,m]; m1[is.na(m1)] <- 0; if(is.null(nrow(m1))) { m1 <- matrix(1,ncol=1) } else {diag(m1)[is.na(m)] <- 1 }
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

saveRDS(summ_list,paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/summ_list_chr",chr,".rds"))
saveRDS(snps_scale,paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/snps_scale_chr",chr,".rds"))
saveRDS(LD_list,paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/LD_list_chr",chr,".rds"))

save(indx, indx_block, delta_best, lambda_best, summ_max, M, snp_list, Nsnps, ethnic, N,
     file=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/chr",chr,".RData"))

