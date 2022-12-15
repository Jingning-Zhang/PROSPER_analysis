rm(list = ls())

args <- commandArgs(T)
for(ii in 1:length(args)){ eval(parse(text=args[[ii]])) }
rep <- as.integer(rep)

library(bigreadr)
library(readr)

L=20
S=10

NCORES=12

load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/EUR/rep_",rep,"/best_tuning_rho_",rho,"_size_4.RData"))

system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/", ethnic,"/tv_allchr",
                " --score /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/EUR/rep_",rep,"/beta_rho_",rho,"_size_4.txt",
                " cols=+scoresums,-scoreavgs --score-col-nums ",best+2,
                " --out /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_4_EUR_best"))
SCORE_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_4_EUR_best.sscore"))
pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",ethnic, "/phenotypes_rho",rho,"_",setting,".phen"))

################################################
################################################

## validation


validation_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/validation.id.txt"), stringsAsFactors = F)$V1
SCORE <- SCORE_all[match(validation_id, SCORE_all$IID), -1:-4]

pheno <- pheno_all[match(validation_id, pheno_all$V2), ]
pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep]

dat <- data.frame(y=pheno, score=SCORE )
fit <- lm(y~., data = dat)
R2 <- summary(fit)$r.squared

saveRDS(R2, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_4_EUR_best.rds"))

R2

