# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################
# closest3lambda_2ethnicity

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
sourceCpp('MEPRS12.2.cpp')


### Load standard data
summ_list <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/summ_list_fixEUR_rho_",rho,"_size_",size,"_chr",chr,".rds"))
LD_list <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/LD_list_rho_",rho,"_size_",size,"_chr",chr,".rds"))
# indx, indx_block, M, snp_list, Nsnps, ethnic, delta_best, summ_max,
load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/otherinfo_rho_",rho,"_size_",size,"_chr",chr,".RData"))

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

m <- which(ethnic=="EUR")

start_time <- Sys.time()
res <- enet_multiethnic(summ=summ_list, R=LD_list,
                        M=M, m = m-1,
                        indx=indx,
                        indx_block=indx_block,
                        delta=delta_best, lambdapath=lambdapath, cpath=cpath)
end_time <- Sys.time()

runtime <- as.numeric(end_time-start_time, units="secs")

save(ethnic, M, delta_best, runtime, res,
     file= paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/4_multi-lassosum_by_chr_",para,"_fixEUR/tuning/rho_",rho,"_size_",size,"_chr",chr,".RData"))

