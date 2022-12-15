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

df_beta <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic,"/rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting,"_matched_to_ref.rds"))

corr <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr/",ethnic,"/corr_sfbm.rds"))
map <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr/",ethnic,"/map.rds"))
df_beta <- df_beta[match(map$id,df_beta$rsid),]
m <- df_beta$a1 != map$a1
if(sum(m)>0){
  df_beta$a1[m] <- map$a1[m]
  df_beta$a0[m] <- map$a0[m]
  df_beta$beta[m] <- -df_beta$beta[m]
}

# ------------------ Run lassosum2 ------------------
NCORES=1
start_time <- Sys.time()
beta_lassosum2 <- snp_lassosum2(corr, df_beta, delta = delta_path(max=100,min=0.01,n=10),
                                ncores = NCORES, maxiter=1000)
params2 <- attr(beta_lassosum2, "grid_param")
beta_lassosum2[is.na(beta_lassosum2)] <- 0
end_time <- Sys.time()
runtime <- as.numeric(end_time-start_time, units="secs")

## the output beta effect sizes are in original scale!!!
save(df_beta, beta_lassosum2, params2, runtime,
#save(df_beta, beta_lassosum2, params2,
     file = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/beta_rho_",rho,"_size_",size,".RData"))
Beta_all <- cbind(df_beta[,c(5,4)], beta_lassosum2)
fwrite2(Beta_all, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/beta_rho_",rho,"_size_",size,".txt"), col.names = F, sep="\t", nThread=NCORES)
write_tsv(params2, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/param_rho_",rho,"_size_",size,".txt"))

