rm(list = ls())

args <- commandArgs(T)
for(ii in 1:length(args)){ eval(parse(text=args[[ii]])) }


library(bigparallelr)
library(bigsnpr)
library(bigreadr)
library(readr)

delta_path <- function (max=100, min=0.5, n=10){
   sqrt_max <- max^(1/3)
   sqrt_min <- min^(1/3)
  path <- numeric(n)
  for (i in 1:n) {
    path[n+1-i] = (sqrt_max - (sqrt_max-sqrt_min)/(n-1)  * (i-1) )^3;
  }
  return(path)
}

#ethnic="EUR"
#trait="depression"

load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/",ethnic,"_",trait,"_corr.RData"))

# ------------------ Run lassosum2 ------------------
NCORES=2
beta_lassosum2 <- snp_lassosum2(corr, df_beta, delta = delta_path(max=100,min=0.5,n=10),
                                ncores = NCORES, maxiter=1000)
params2 <- attr(beta_lassosum2, "grid_param")
beta_lassosum2[is.na(beta_lassosum2)] <- 0

save(df_beta, beta_lassosum2, params2, file = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/",ethnic,"_",trait,".RData"))

