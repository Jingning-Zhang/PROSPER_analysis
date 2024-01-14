
rm(list=ls())

library(bigreadr)
library(readr)
library(dplyr)


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

NCORES <-  1

setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",race,"/rep_",rep))

dir.create(paste0("./ldpred2/prs/results/"), recursive=T)


################################################
################################################
## validation data

pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",race, "/phenotypes_rho",rho,"_",setting,".phen"))

validation_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/validation.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(validation_id, pheno_all$V2), ]

race1 <- "EUR"
best <- readRDS(paste0("../../",race1,"/rep_",rep,"/ldpred2/prs/results/rho_",rho,"_size_",ifelse(race1=='EUR',4,size),"/",race1,"_best_param.rds"))

load(paste0("./ldpred2/prs/validation/rho_",rho,"_size_",size,"/",race,"-",race1,"_validation_score.RData"))
SCORE <- SCORE[,best,drop=F]

iid <- intersect(pheno$V2, SCORE_id$IID)

pheno <- pheno[match(iid, pheno$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ,drop=F]

pheno_validation <- pheno
SCORE_validation <- SCORE

################################################
## test

fit = lm(pheno_validation~SCORE_validation)
R2 <- summary(fit)$r.square

tmp <- data.frame(EUR_R2=R2)
res <- tmp

saveRDS(res, paste0("./ldpred2/prs/results/rho_",rho,"_size_",size,"/", race,"_validation_R2_EUR_best.rds"))
res



