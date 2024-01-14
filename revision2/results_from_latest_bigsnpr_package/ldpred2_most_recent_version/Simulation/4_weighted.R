# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())


library(bigreadr)
library(readr)
library(dplyr)
#library(SuperLearner)

parameters <- expand.grid(race = c('AFR','AMR','EAS','SAS'),
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

NCORES = 1


################################################
## prepare data
################################################

## score

race1_all <- c('AFR','AMR','EAS','EUR','SAS')

pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",race, "/phenotypes_rho",rho,"_",setting,".phen"))

## best tuning parameters
best <- NULL
for(i in 1:length(race1_all) ){
  race1 <- race1_all[i]
  best <- c(best, readRDS(paste0("../../",race1,"/rep_",rep,"/ldpred2/prs/results/rho_",rho,"_size_",ifelse(race1=='EUR',4,size),"/",race1,"_best_param.rds")))
}
names(best) <- race1_all

################################################
## tuning data

tuning_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/test.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(tuning_id, pheno_all$V2), ]

for(i in 1:length(race1_all) ){
  race1 <- race1_all[i]
  load(paste0("./ldpred2/prs/tuning/rho_",rho,"_size_",size,"/",race,"-",race1,"_tuning_score.RData"))
  if(race1=='AFR'){
    id <- SCORE_id
    SCORE_all <- SCORE
    SCORE_best <- SCORE[,best[race1]]
  }else{
    SCORE_all <- cbind(SCORE_all, SCORE[match(id$IID, SCORE_id$IID),])
    SCORE_best <- cbind(SCORE_best, SCORE[match(id$IID, SCORE_id$IID),best[race1]])
  }
}
SCORE_id <- id
SCORE <- SCORE_all

iid <- intersect(pheno$V2, SCORE_id$IID)

pheno <- pheno[match(iid, pheno$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
SCORE_best <- as.matrix(SCORE_best)[match(iid, SCORE_id$IID), ]

pheno_tuning <- pheno
SCORE_tuning <- SCORE
SCORE_tuning_best <- SCORE_best

################################################
## validation data

validation_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/validation.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(validation_id, pheno_all$V2), ]

for(i in 1:length(race1_all) ){
  race1 <- race1_all[i]
  load(paste0("./ldpred2/prs/validation/rho_",rho,"_size_",size,"/",race,"-",race1,"_validation_score.RData"))
  if(race1=='AFR'){
    id <- SCORE_id
    SCORE_all <- SCORE
    SCORE_best <- SCORE[,best[race1]]
  }else{
    SCORE_all <- cbind(SCORE_all, SCORE[match(id$IID, SCORE_id$IID),])
    SCORE_best <- cbind(SCORE_best, SCORE[match(id$IID, SCORE_id$IID),best[race1]])
  }
}
SCORE_id <- id
SCORE <- SCORE_all

iid <- intersect(pheno$V2, SCORE_id$IID)

pheno <- pheno[match(iid, pheno$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
SCORE_best <- as.matrix(SCORE_best)[match(iid, SCORE_id$IID), ]

pheno_validation <- pheno
SCORE_validation <- SCORE
SCORE_validation_best <- SCORE_best

################################################
## test

## ridge
#
#tmp <- data.frame(SCORE_tuning); m <- apply(tmp,2,sd)==0; tmp <- tmp[,!m, drop=F]
#sl = SuperLearner(Y = pheno_tuning,
#                  X = tmp,
#                  family = gaussian(),
#                  SL.library = c("SL.ridge") )
#
#score <- predict(sl, data.frame(SCORE_validation)[,!m, drop=F], onlySL = TRUE)[[1]]
#
#fit = lm(pheno_validation~score)
#R2 <- summary(fit)$r.square
#tmp <- data.frame(ridge_combined_R2=R2)
#res <- tmp
#sl1 = sl
#
#
## ridge+lasso+lm
#tmp <- data.frame(SCORE_tuning); m <- apply(tmp,2,sd)==0; tmp <- tmp[,!m, drop=F]
#sl = SuperLearner(Y = pheno_tuning,
#                  X = tmp,
#                  family = gaussian(),
#                  SL.library = c("SL.glmnet", "SL.ridge","SL.lm") )
#
#score <- predict(sl, data.frame(SCORE_validation)[,!m, drop=F], onlySL = TRUE)[[1]]
#
#fit = lm(pheno_validation~score)
#R2 <- summary(fit)$r.square
#tmp <- data.frame(superlearning_combined_R2=R2)
#res <- cbind(res,tmp)
#
#sl2 = sl

# weighted

tmp <- data.frame(pheno_tuning,SCORE_tuning_best); colnames(tmp) <- c("y",race1_all)
wt = lm(y~., tmp)
tmp <- data.frame(SCORE_validation_best); colnames(tmp) <- race1_all
score <- predict(wt, tmp)

fit = lm(pheno_validation~score)
R2 <- summary(fit)$r.square
tmp <- data.frame(weighted_combined_R2=R2)
#res <- cbind(res,tmp)
res <- tmp

saveRDS(res, paste0("./ldpred2/prs/results/rho_",rho,"_size_",size,"/",race,"_validation_R2_sl.rds"))

#res_model <- list(ridge=sl1, sl=sl2, wt = wt)
#saveRDS(res_model, paste0("./ldpred2/prs/results/rho_",rho,"_size_",size,"/",race,"_sl_models.rds"))

res

