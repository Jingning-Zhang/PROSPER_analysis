
rm(list=ls())

library(bigsnpr)
library(bigreadr)
library(bigparallelr)
library(stringr)
library(dplyr)


parameters <- expand.grid(race = c('AFR','AMR','EAS','EUR','SAS'),
                          setting = 1:5, rho = 1:3, size = 1:4, rep = 1:10)


args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = as.character(parameters[as.numeric(temp1),1])
setting = as.integer(parameters[as.numeric(temp1),2])
rho = as.integer(parameters[as.numeric(temp1),3])
size = as.integer(parameters[as.numeric(temp1),4])
rep = as.integer(parameters[as.numeric(temp1),5])

print(setting)
print(rho)
print(size)
print(rep)

print(race)

NCORES <- 1

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",race,"/rep_",rep), recursive=T)
setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",race,"/rep_",rep))

dir.create("./lassosum2")
dir.create("./lassosum2/coef")
dir.create("./lassosum2/prs")
dir.create("./lassosum2/prs/tuning")

df_beta <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",race,"/rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting,"_matched_to_ref.rds"))

corr <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr/",race,"/corr_sfbm.rds"))
map <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr/",race,"/map.rds"))
df_beta <- df_beta[match(map$id,df_beta$rsid),]
m <- df_beta$a1 != map$a1
if(sum(m)>0){
  df_beta$a1[m] <- map$a1[m]
  df_beta$a0[m] <- map$a0[m]
  df_beta$beta[m] <- -df_beta$beta[m]
}

set.seed(1)  # to get the same result every time
beta_grid <- snp_lassosum2(corr, df_beta, ncores = NCORES)
save(beta_grid, df_beta, file=paste0("./lassosum2/coef/rho_",rho,"_size_",size,"_lassosum2.RData"))

load(paste0("./lassosum2/coef/rho_",rho,"_size_",size,"_lassosum2.RData"))

beta_grid[is.na(beta_grid)] = 0
beta_grid = data.frame(df_beta[,c('rsid','a1')], beta_grid)
colnames(beta_grid) = c('rsid','a1', paste0('e',1:(ncol(beta_grid)-2)))

for(chr in 1:22){
  print(chr)

  m <- df_beta$chr==chr
  prs.file <- beta_grid[m,,drop=F]
  write.table(prs.file, file = paste0("./lassosum2/coef/rho_",rho,"_size_",size,"_lassosum2-chr",chr,".txt"),
              col.names = F,row.names = F,quote=F)

}


