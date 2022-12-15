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

## prs calculation

#for (mm in 1:M){
#
#  for (chr in 1:22){
#
#    print(paste0(mm,"-",chr))
#
#    system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
#                  " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/",targetethnic,"/chr", chr,
#                  " --score ", path,"/Results/ukbb_beta_",ethnic[mm],"_",trait,".txt",
#                  " cols=+scoresums,-scoreavgs --score-col-nums 3-",alltuning[[mm]]+2,
#                  " --out ",path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_",trait,"_chr", chr,"_tuning_score"))
#
#    score <- fread2(paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_",trait,"_chr", chr,"_tuning_score.sscore"))
#    if(chr==1){
#      SCORE_all <- score[, -1:-4,drop=F]
#      SCORE_id <- score[,1:2,drop=F]
#    }else{
#      SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
#    }
#  }
#  SCORE <- SCORE_all
#  save(SCORE, SCORE_id, file = paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_",trait,"_tuning_score.RData"))
#
#}
#
#for (mm in 1:M){
#  load(paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_",trait,"_tuning_score.RData"))
#  if(mm==1){
#    SCORE_all <- SCORE
#    SCORE_id_all <- SCORE_id
#  }else{
#    SCORE_all <- cbind(SCORE_all, SCORE[match(SCORE_id_all$IID, SCORE_id$IID),])
#  }
#}
#SCORE_id <- SCORE_id_all
#SCORE <- SCORE_all
#save(SCORE, SCORE_id, prseth,
#     file = paste0(path,"/Results/test_all_prs/",targetethnic,"/all_ancestry_",trait,"_tuning_score.RData"))

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

## prs calculation

#for (mm in 1:M){
#
#  for (chr in 1:22){
#
#    print(paste0(mm,"-",chr))
#
#    system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
#                  " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/",targetethnic,"/chr", chr,
#                  " --score ", path,"/Results/ukbb_beta_",ethnic[mm],"_",trait,".txt",
#                  " cols=+scoresums,-scoreavgs --score-col-nums 3-",alltuning[[mm]]+2,
#                  " --out ",path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_",trait,"_chr", chr,"_validation_score"))
#
#    score <- fread2(paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_",trait,"_chr", chr,"_validation_score.sscore"))
#    if(chr==1){
#      SCORE_all <- score[, -1:-4,drop=F]
#      SCORE_id <- score[,1:2,drop=F]
#    }else{
#      SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
#    }
#  }
#  SCORE <- SCORE_all
#  save(SCORE, SCORE_id, file = paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_",trait,"_validation_score.RData"))
#
#}
#
#for (mm in 1:M){
#  load(paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_",trait,"_validation_score.RData"))
#  if(mm==1){
#    SCORE_all <- SCORE
#    SCORE_id_all <- SCORE_id
#  }else{
#    SCORE_all <- cbind(SCORE_all, SCORE[match(SCORE_id_all$IID, SCORE_id$IID),])
#  }
#}
#SCORE_id <- SCORE_id_all
#SCORE <- SCORE_all
#save(SCORE, SCORE_id, prseth,
#     file = paste0(path,"/Results/test_all_prs/",targetethnic,"/all_ancestry_",trait,"_validation_score.RData"))

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

## prs evaluation

#############################################################
#############################################################

#################################
# best tuning prs

## tuning

R2 <- numeric(length = nprs)
for (i in 1:nprs){
  model = lm(y~.-id,data = pheno_tuning); residual = model$residuals
  fit <- lm( residual~SCORE_tuning[,i] )
  R2[i] <- summary(fit)$r.square
}
best <- which.max(R2)
best_target <- which(prseth==targetethnic)[which.max(R2[prseth==targetethnic])]

tmp <- data.frame(best=best_target, best_across_ancestry=best)

saveRDS(tmp, paste0(path,"/Results/final_results/",targetethnic,"/",trait,"_best.rds"))

## validation

score <- SCORE_validation[,best]
model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit <- lm( residual~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(best_R2=R2)
res <- tmp

score <- SCORE_validation[,best_target]
model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit <- lm( residual~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(best_R2_across_ancestry=R2)
res <- cbind(res,tmp)

res


#################################
# combined prs

## lasso-combined target ancestry

set.seed(20220716)

model = lm(y~.-id,data = pheno_tuning)
residual = model$residuals

a <- SCORE_tuning[,prseth==targetethnic,drop=F]
b <- SCORE_validation[,prseth==targetethnic,drop=F]
m <- apply(a,2,sd)
m <- which(is.na(m) | m==0)
if(length(m)>0){
  a <- a[,-m,drop=F]
  b <- b[,-m,drop=F]
}
  sl = SuperLearner(Y = residual,
                    X = data.frame(a),
                    family = gaussian(),
                    SL.library = c("SL.glmnet") )

score <- predict(sl, data.frame(b), onlySL = TRUE)[[1]]

model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit = lm(residual~score)
R2 <- summary(fit)$r.square

tmp <- data.frame(lasso_combined_R2=R2)
res <- cbind(res,tmp)

sl1=sl


## lasso-combined all ancestry

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
                    SL.library = c("SL.glmnet") )

score <- predict(sl, data.frame(b), onlySL = TRUE)[[1]]

model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit = lm(residual~score)
R2 <- summary(fit)$r.square

tmp <- data.frame(lasso_combined_R2_cross_ancestry=R2)
res <- cbind(res,tmp)

sl2=sl

saveRDS(res, paste0(path,"/Results/final_results/",targetethnic,"/",trait,"_validation_R2_no_sl.rds"))

res


## superlearning-combined target ancestry

set.seed(20220716)

model = lm(y~.-id,data = pheno_tuning)
residual = model$residuals

a <- SCORE_tuning[,prseth==targetethnic,drop=F]
b <- SCORE_validation[,prseth==targetethnic,drop=F]
m <- apply(a,2,sd)
m <- which(is.na(m) | m==0)
if(length(m)>0){
  a <- a[,-m,drop=F]
  b <- b[,-m,drop=F]
}

  sl = SuperLearner(Y = residual,
                    X = data.frame(a),
                    family = gaussian(),
                    SL.library = c("SL.glmnet", "SL.ridge", "SL.lm") )

score <- predict(sl, data.frame(b), onlySL = TRUE)[[1]]

model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit = lm(residual~score)
R2 <- summary(fit)$r.square

tmp <- data.frame(superlearning_combined_R2=R2)
res <- cbind(res,tmp)

sl3=sl

## superlearning-combined all ancestry

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
                    SL.library = c("SL.glmnet", "SL.ridge", "SL.lm") )

score <- predict(sl, data.frame(b), onlySL = TRUE)[[1]]

model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit = lm(residual~score)
R2 <- summary(fit)$r.square

tmp <- data.frame(superlearning_combined_R2_cross_ancestry=R2)
res <- cbind(res,tmp)


sl4=sl

saveRDS(res, paste0(path,"/Results/final_results/",targetethnic,"/",trait,"_validation_R2.rds"))


res_sl <- list(lasso=sl1, lasso_across_ancestry=sl2, sl=sl3, sl_across_ancestry=sl4)


saveRDS(res_sl, paste0(path,"/Results/final_results/",targetethnic,"/",trait,"_combined_models.rds"))

res

