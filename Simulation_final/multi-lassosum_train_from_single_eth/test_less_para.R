# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

#setting=2
#rep=1
#rho=1
#size=1
#
#ethnic1="AFR_EUR"
#targetethnic="AFR"
#para=10

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
ethnic1 <- paste(ethnic,collapse = "_")


progBar <- function(ii, N, per = 10) {
  #ii is current iteration.
  #N is total number of iterations to perform.
  #per is step (percent) when to update the progress bar. We need this as when multiple iterations are being performed progress bar gets updated too often and output is messed up.
  if (ii %in% seq(1, N, per)) {
    x <- round(ii * 100 / N)
    message("[ ",
            paste(rep("=", x), collapse = ""),
            paste(rep("-", 100 - x), collapse = ""),
            " ] ", x, "%", "\r",
            appendLF = FALSE)
    if (ii == N) cat("\r")
  }
}

library(doMC)
library(foreach)
registerDoMC(6)

library(bigreadr)
library(readr)
library(dplyr)
library(glmnet)
#library(caret)

if(fixEUR=="T"){
  path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/4_multi-lassosum_by_chr_",para,"_fixEUR")
}else{
  path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para)
}


pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",targetethnic, "/phenotypes_rho",rho,"_",setting,".phen"))

load(paste0(path, "/Results/",targetethnic,"_rho_",rho,"_size_",size,"_tv_score.RData"))

parameter <- read_tsv(paste0(path, "/Results/para_rho_",rho,"_size_",size,".txt"))

#lam1 <- unique(parameter$lambda1)[c(1,3)]
#lam2 <- unique(parameter$lambda2)[c(1,3)]
#lam3 <- unique(parameter$lambda3)[c(1,3)]
#lam4 <- unique(parameter$lambda4)[c(1,3)]
#lam5 <- unique(parameter$lambda5)[c(1,3)]
#c <- unique(parameter$c)[c(1,3)]
#mm <- which( (parameter$lambda1 %in% lam1) & (parameter$lambda2 %in% lam2) & (parameter$lambda3 %in% lam3) & (parameter$lambda4 %in% lam4) & (parameter$lambda5 %in% lam5) & (parameter$c %in% c) )

lam1 <- unique(parameter$lambda1)[c(5,10)]
lam2 <- unique(parameter$lambda2)[c(5,10)]
c <- unique(parameter$c)[c(1,10)]
mm <- which( (parameter$lambda1 %in% lam1) & (parameter$lambda2 %in% lam2) & (parameter$c %in% c) )

## tuning
id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/geno/mega_matched_summary_stat/",targetethnic,"/test.id.txt"), stringsAsFactors = F)$V1

pheno <- pheno_all[match(id, pheno_all$V2), ]; pheno <- as.matrix(pheno[,-1:-2])
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,] }

#mtx = cor(SCORE); drop = findCorrelation(mtx,cutoff = 0.98)
#SCORE_tuning = SCORE[,-drop]
SCORE_tuning = SCORE[,mm]
pheno_tuning <- pheno
#rm(list = c("mtx","pheno","SCORE"))



## validation
id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_multiethnic-lassosum2-parameter/geno/mega_matched_summary_stat/",targetethnic,"/validation.id.txt"), stringsAsFactors = F)$V1

pheno <- pheno_all[match(id, pheno_all$V2), ]; pheno <- as.matrix(pheno[,-1:-2])
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,] }
pheno_validation <- pheno
#SCORE_validation = SCORE[,-drop]
SCORE_validation = SCORE[,mm]

rm(list = c("pheno","SCORE"))

## lassoSCORE
set.seed(20220210)
niter=100
R2 <- foreach(i=1:niter, ii = icount(), .combine='c') %dopar% {
  fit <- glmnet::cv.glmnet(x=as.matrix(SCORE_tuning), y=pheno_tuning[,i], family='gaussian', alpha=1)
  score <- as.numeric(predict(fit, as.matrix(SCORE_validation), s='lambda.min'))
  fit <- lm(pheno_validation[,i] ~ score)
  progBar(ii, niter, per=5)
  return(summary(fit)$r.squared)
}

tmp <- is.null(R2)
if(sum(tmp) != 0){
  for (i in which(tmp) ){
    fit <- glmnet::cv.glmnet(x=as.matrix(SCORE_tuning), y=pheno_tuning[,i], family='gaussian', alpha=1)
    score <- as.numeric(predict(fit, as.matrix(SCORE_validation), s='lambda.min'))
    fit <- lm(pheno_validation[,i] ~ score)
    R[i] <- summary(fit)$r.squared
  }
}

mean(R2)

#mean(readRDS(paste0(path,"/Results/",targetethnic,"_rho_",rho,"_size_",size,"_validation_lassoSCORE_R2_all.rds")))
mean(readRDS(paste0(path,"/Results/",targetethnic,"_rho_",rho,"_size_",size,"_validation_superlearning_R2_all.rds")))
mean(readRDS(paste0(path,"/Results/",targetethnic,"_rho_",rho,"_size_",size,"_validation_R2.rds")))


