
rm(list=ls())

library(bigsnpr)
library(bigreadr)
library(bigparallelr)
library(readr)
library(dplyr)
library(SuperLearner)


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


setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",race,"/rep_",rep))

dir.create(paste0("./lassosum2/prs/results/rho_",rho,"_size_",size), recursive=T)

pheno_all = read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",race,"/phenotypes_rho",rho,"_",setting,".phen"))

################################################
################################################

## tuning

load(paste0("./lassosum2/prs/tuning/rho_",rho,"_size_",size,"/",race,"-",race,"_tuning_score.RData"))

tuning_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/test.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(tuning_id, pheno_all$V2), ]

iid <- intersect(pheno$V2, SCORE_id$IID)

pheno <- pheno[match(iid, pheno$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]

pheno_tuning <- pheno
SCORE_tuning <- SCORE

alltuning <- ncol(SCORE)

R2 <- numeric(length = alltuning)
for (i in 1:(alltuning)){
  fit <- lm( pheno ~ SCORE[,i] )
  R2[i] <- summary(fit)$r.square
}
best <- which.max(R2)
saveRDS(best, paste0("./lassosum2/prs/results/rho_",rho,"_size_",size,"/",race,"_best_param.rds"))


################################################
################################################

## validation

load(paste0("./lassosum2/prs/validation/rho_",rho,"_size_",size,"/",race,"-",race,"_validation_score.RData"))

validation_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/validation.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(validation_id, pheno_all$V2), ]

iid <- intersect(pheno$V2, SCORE_id$IID)

pheno <- pheno[match(iid, pheno$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]

pheno_validation <- pheno
SCORE_validation <- SCORE

score <- SCORE_validation[,best]

fit <- lm( pheno~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(best_R2=R2)
res <- tmp

saveRDS(res, paste0("./lassosum2/prs/results/rho_",rho,"_size_",size,"/",race,"_validation_R2_best.rds"))
res
