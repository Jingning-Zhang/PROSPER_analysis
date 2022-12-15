# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)


library(dplyr)
library(glmnet)

library(SuperLearner)
library(caret)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
L <- as.numeric(L); Lc <- as.numeric(Lc)

NCORES <- 1

path="/dcs04/nilanjan/data/23andme/Analysis/JZ/weighted_lassosum2/ukbb_test/"

#############################################################
#############################################################

## two ethnic

  for (chr in 1:22){

    print(paste0(chr))

    system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                  " --bfile /dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/",ethnic,"/chr", chr,
                  " --score ",path,"/",ethnic,"_",trait,".prs",
                  #" cols=+scoresums,-scoreavgs --score-col-nums 3", ifelse( (fixEUR == "T") & (ethnic[mm]=="EUR"), "", paste0("-",alltuning+2)),
                  " cols=+scoresums,-scoreavgs --score-col-nums 3-4",
                  " --out ",path,"/",ethnic,"_",trait,"_chr", chr,"_tv_score"))

    score <- fread2(paste0(path,"/",ethnic,"_",trait,"_chr", chr,"_tv_score.sscore"))
    if(chr==1){
      SCORE_all <- score[, -1:-4,drop=F]
      SCORE_id <- score[,1:2,drop=F]
    }else{
      SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
    }
  }
  SCORE <- SCORE_all
  save(SCORE, SCORE_id, file = paste0(path,"/",ethnic,"_",trait,"_tv_score.RData"))


###########
pheno_all = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/',ethnic,'/pheno.txt'),header=T)
pheno_all = pheno_all[,c(trait,'eid','age','sex',paste0('pc',1:10))]
colnames(pheno_all) = c('y','id','age','sex',paste0('pc',1:10))
pheno_all = pheno_all[complete.cases(pheno_all$y),]
pheno_all$id = as.character(pheno_all$id)


## tuning

id = read.table(paste0('/dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/',ethnic,'/bmi.IDs.tuning.txt'),header=F)$V2
id <- intersect(id, intersect(pheno_all$id, SCORE_id$IID))
pheno <- pheno_all[match(id, pheno_all$id), ]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,] }
SCORE_tuning = SCORE; colnames(SCORE_tuning) = c("EUR","TAR")
pheno_tuning <- pheno


## validation

id = read.table(paste0('/dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/',ethnic,'/bmi.IDs.validation.txt'),header=F)$V2
pheno <- pheno_all[match(id, pheno_all$id), ]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,] }
SCORE_validation = SCORE; colnames(SCORE_validation) = c("EUR","TAR")
pheno_validation <- pheno
rm(list = c("pheno","SCORE"))



set.seed(20220502)
model = lm(y~.-id,data = pheno_tuning)
residual = model$residuals

model2 = lm(residual~SCORE_tuning)
weight = coef(model2)[2:3]
PRS = SCORE_validation%*%weight

model1 = lm(y~.-id,data = pheno_validation)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square

# 0.06618649



#############################################################
#############################################################

## all ethnic

  for (chr in 1:22){

    print(paste0(chr))

    system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                  " --bfile /dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/",ethnic,"/chr", chr,
                  " --score ",path,"/all_ethnic_",trait,".prs",
                  #" cols=+scoresums,-scoreavgs --score-col-nums 3", ifelse( (fixEUR == "T") & (ethnic[mm]=="EUR"), "", paste0("-",alltuning+2)),
                  " cols=+scoresums,-scoreavgs --score-col-nums 3-7",
                  " --out ",path,"/all_ethnic_",ethnic,"_",trait,"_chr", chr,"_tv_score"))

    score <- fread2(paste0(path,"/all_ethnic_",ethnic,"_",trait,"_chr", chr,"_tv_score.sscore"))
    if(chr==1){
      SCORE_all <- score[, -1:-4,drop=F]
      SCORE_id <- score[,1:2,drop=F]
    }else{
      SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
    }
  }
  SCORE <- SCORE_all
  save(SCORE, SCORE_id, file = paste0(path,"/all_ethnic_",ethnic,"_",trait,"_tv_score.RData"))


################################################
pheno_all = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/',ethnic,'/pheno.txt'),header=T)
pheno_all = pheno_all[,c(trait,'eid','age','sex',paste0('pc',1:10))]
colnames(pheno_all) = c('y','id','age','sex',paste0('pc',1:10))
pheno_all = pheno_all[complete.cases(pheno_all$y),]
pheno_all$id = as.character(pheno_all$id)

################################################
################################################

## tuning

id = read.table(paste0('/dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/',ethnic,'/bmi.IDs.tuning.txt'),header=F)$V2
id <- intersect(id, intersect(pheno_all$id, SCORE_id$IID))
pheno <- pheno_all[match(id, pheno_all$id), ]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,] }
SCORE_tuning = SCORE; colnames(SCORE_tuning) = c("EUR","AFR","AMR","EAS","SAS")
pheno_tuning <- pheno


## validation

id = read.table(paste0('/dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/',ethnic,'/bmi.IDs.validation.txt'),header=F)$V2
pheno <- pheno_all[match(id, pheno_all$id), ]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,] }
SCORE_validation = SCORE; colnames(SCORE_validation) = c("EUR","AFR","AMR","EAS","SAS")
pheno_validation <- pheno
rm(list = c("pheno","SCORE"))



set.seed(20220502)
model = lm(y~.-id,data = pheno_tuning)
residual = model$residuals

model2 = lm(residual~SCORE_tuning)
weight = coef(model2)[2:6]
PRS = SCORE_validation%*%weight

model1 = lm(y~.-id,data = pheno_validation)
residual = model1$residuals
model2 = lm(residual~PRS)
summary(model2)$r.square

# 0.07504062

