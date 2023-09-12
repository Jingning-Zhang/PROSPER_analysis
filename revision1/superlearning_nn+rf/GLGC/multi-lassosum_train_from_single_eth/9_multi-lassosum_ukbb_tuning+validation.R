## /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
###############################
#
#rm(list=ls())

library(readr)
library(bigreadr)

library(dplyr)
library(glmnet)

library(SuperLearner)
library(caret)

#args <- commandArgs(T)
#for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
L <- as.numeric(L); Lc <- as.numeric(Lc)

NCORES <- 1

#############################################################
#############################################################

M <- length(ethnic)
ethnic1 <- paste(ethnic,collapse = "_")

path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para)


alltuning <- list()
for (mm in 1:M){
  alltuning[[mm]] <- readRDS(paste0(path,"/Results/alltuning_",ethnic[mm],"_",trait,".rds"))
}

#############################################################
#############################################################


nprs <- sum(unlist(alltuning))
prseth  <- character()
for (mm in 1:M){
  prseth <- c(prseth, rep(ethnic[mm], alltuning[[mm]]))
}


#############################################################
#############################################################

## tuning

#############################################################
#############################################################

load(paste0(path,"/Results/test_all_prs/",targetethnic,"/all_ancestry_",trait,"_tuning_score.RData")) #!!!

## phenotype

pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",targetethnic,"_tuning.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",targetethnic,"_tuning.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

iid <- intersect(pheno$id, SCORE_id$IID)
pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
pheno_tuning <- pheno
SCORE_tuning <- SCORE

#############################################################
#############################################################

## validation

#############################################################
#############################################################

load(paste0(path,"/Results/test_all_prs/",targetethnic,"/all_ancestry_",trait,"_validation_score.RData")) #!!!

## phenotype

pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",targetethnic,"_validation.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",targetethnic,"_validation.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

iid <- intersect(pheno$id, SCORE_id$IID)
pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
pheno_validation <- pheno
SCORE_validation <- SCORE

#############################################################
#############################################################

## combine scores using other machine learning methods

#############################################################
#############################################################


set.seed(20220716)

model = lm(y~.-id,data = pheno_tuning)
residual = model$residuals

a <- SCORE_tuning
b <- SCORE_validation
m <- apply(a,2,sd)
m <- which(is.na(m) | m==0)
if(length(m)>0){
  a <- a[,-m,drop=F]
  b <- b[,-m,drop=F]
}

  sl = SuperLearner(Y = residual,
                    X = data.frame(a),
                    family = gaussian(),
                    SL.library = c("SL.glmnet", "SL.ridge", "SL.lm","SL.nnet", "SL.randomForest") )

score <- predict(sl, data.frame(b), onlySL = TRUE)[[1]]

model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit = lm(residual~score)
R2 <- summary(fit)$r.square

tmp <- data.frame(superlearning_combined_R2_cross_ancestry=R2)

tmp
