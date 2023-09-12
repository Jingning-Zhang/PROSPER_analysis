# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

library(doMC)
library(foreach)


library(dplyr)
library(glmnet)

library(SuperLearner)
library(caret)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
L <- as.numeric(L); Lc <- as.numeric(Lc)
rep <- as.integer(rep)
Ntuning <- as.integer(Ntuning)

NCORES <- 1
registerDoMC(NCORES)

#############################################################
#############################################################

M <- length(ethnic)
ethnic1 <- paste(ethnic,collapse = "_")


path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/number_of_tuning_samples/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para)

whichethnic <- which(ethnic == targetethnic)

load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para,"/Results/test_all_prs/",targetethnic,"/all_prs_rho_",rho,"_size_",size,"_tv_score.RData"))

pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",targetethnic, "/phenotypes_rho",rho,"_",setting,".phen"))


################################################
################################################

## combined prs

################################################
################################################

## tuning

id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",targetethnic,"/test.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(id, pheno_all$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,,drop=F] }
pheno_tuning <- pheno
SCORE_tuning = SCORE


## validation
id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",targetethnic,"/validation.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(id, pheno_all$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,,drop=F] }
pheno_validation <- pheno
SCORE_validation = SCORE
rm(list = c("pheno","SCORE"))


set.seed(20230823)
mm <- sort(sample(1:length(pheno_tuning), Ntuning))

set.seed(20220502)

a <- SCORE_tuning
b <- SCORE_validation
m <- apply(a,2,sd)
m <- which(is.na(m) | m==0)
if(length(m)>0){
  a <- a[,-m,drop=F]
  b <- b[,-m,drop=F]
}

  sl = SuperLearner(Y = as.numeric(pheno_tuning[mm]),
                    X = data.frame(a[mm,]),
                    family = gaussian(),
                    SL.library = c("SL.glmnet", "SL.ridge","SL.lm") )
  score = predict(sl, data.frame(b), onlySL = TRUE)[[1]]
  fit <- lm(as.numeric(pheno_validation) ~ score)
  R2 <- summary(fit)$r.square

res <- data.frame(superlearning_combined_R2_cross_ancestry=R2)

saveRDS(res, paste0(path,"/Results/test_all_prs/",targetethnic,"/rho_",rho,"_size_",size,"_Ntuning",Ntuning,"_validation_R2.rds"))

res


