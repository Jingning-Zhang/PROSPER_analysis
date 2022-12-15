# /dcs04/nilanjan/data/jzhang2/MEPRS
##############################

rm(list=ls())

library(readr)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(inline)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
L <- as.integer(L); Lc <- as.integer(Lc);

M <- length(ethnic)
ethnic1 <- paste(ethnic,collapse = "_")

setwd("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/")
sourceCpp('MEPRS12.1.cpp')

### Load standard data
# indx, indx_block, delta_best, delta_best, summ_max, M, snp_list, Nsnps, ethnic, N,
load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/chr",chr,".RData"))

summ_list <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/summ_list_chr",chr,".rds"))
snps_scale <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/snps_scale_chr",chr,".rds"))
LD_list <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum/",ethnic1,"/1_standard_data/",trait,"/LD_list_chr",chr,".rds"))

##############################
### tuning parameters

lambdapath <- matrix(nrow = M, ncol = L)
for (l in 1:M){
  lambdapath[l,] <- r_path(rmax = min(summ_max/lambda_best), rmin = min(0.001/lambda_best), nr=L ) * lambda_best[l]
}

cpath_tmp <- c_path(maxc=100, minc=2, nc=Lc)
cpath <- list()
for (lc in 1:Lc){
  cpath[[lc]] <- matrix(cpath_tmp[lc],ncol = M, nrow = M)
}

##############################
### run algorithm

start_time <- Sys.time()
res <- enet_multiethnic(summ=summ_list, R=LD_list,
                        M=M, indx=indx,
                        indx_block=indx_block,
                        delta=delta_best, lambdapath=lambdapath, cpath=cpath)
end_time <- Sys.time()

runtime <- as.numeric(end_time-start_time, units="secs")

save(ethnic, M, delta_best, runtime, res,
     file= paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/tuning/",trait,"_chr",chr,".RData"))



