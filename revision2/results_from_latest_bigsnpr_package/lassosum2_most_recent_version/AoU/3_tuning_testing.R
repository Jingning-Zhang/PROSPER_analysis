
rm(list=ls())

library(bigsnpr)
library(bigreadr)
library(bigparallelr)
library(readr)
library(dplyr)
library(SuperLearner)

parameters = expand.grid(race = c('AFR','AMR','EUR'),
                         trait = c('height', 'bmi'))

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = as.character(parameters[as.numeric(temp1),1])
trait = as.character(parameters[as.numeric(temp1),2])

print(trait)
print(race)

#NCORES <-  3

setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait))

dir.create(paste0("./lassosum2/prs/results/"))

################################################
################################################

## tuning

load(paste0("./lassosum2/prs/tuning/",race,"_tuning_score.RData"))

pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",race,"_tuning.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",race,"_tuning.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

iid <- intersect(pheno$id, SCORE_id$IID)
pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
pheno_tuning <- pheno
SCORE_tuning <- SCORE

alltuning <- ncol(SCORE)

R2 <- numeric(length = alltuning)
for (i in 1:(alltuning)){
  model = lm(y~.-id,data = pheno); residual = model$residuals
  fit <- lm( residual ~ SCORE[,i] )
  R2[i] <- summary(fit)$r.square
}
best <- which.max(R2)
saveRDS(best, paste0("./lassosum2/prs/tuning/",race,"_best.rds"))


################################################
################################################

## validation

load(paste0("./lassosum2/prs/validation/",race,"_validation_score.RData"))

pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",race,"_validation.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",race,"_validation.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

iid <- intersect(pheno$id, SCORE_id$IID)
pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
pheno_validation <- pheno
SCORE_validation <- SCORE

score <- SCORE_validation[,best]

model = lm(y~.-id,data = pheno); residual = model$residuals
fit <- lm( residual~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(best_R2=R2)
res <- tmp

saveRDS(res, paste0("./lassosum2/prs/results/",race,"_validation_R2_best.rds"))
res

################################################
################################################

## combined prs

set.seed(20220712)
model = lm(y~.-id,data = pheno_tuning)
residual = model$residuals

  sl = SuperLearner(Y = residual,
                    X = data.frame(SCORE_tuning),
                    family = gaussian(),
                    SL.library = c("SL.glmnet") )

score <- predict(sl, data.frame(SCORE_validation), onlySL = TRUE)[[1]]

model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit = lm(residual~score)
R2 <- summary(fit)$r.square

tmp <- data.frame(lasso_combined_R2=R2)
res <- cbind(res,tmp)

sl1 = sl

set.seed(20220712)
model = lm(y~.-id,data = pheno_tuning)
residual = model$residuals

  sl = SuperLearner(Y = residual,
                    X = data.frame(SCORE_tuning),
                    family = gaussian(),
                    SL.library = c("SL.glmnet", "SL.ridge", "SL.nnet","SL.lm") )

score <- predict(sl, data.frame(SCORE_validation), onlySL = TRUE)[[1]]

model = lm(y~.-id,data = pheno_validation); residual = model$residuals
fit = lm(residual~score)
R2 <- summary(fit)$r.square
tmp <- data.frame(superlearning_combined_R2=R2)
res <- cbind(res,tmp)

sl2 = sl

saveRDS(res, paste0("./lassosum2/prs/results/",race,"_validation_R2_sl.rds"))

res_sl <- list(lasso=sl1, sl=sl2)
saveRDS(res_sl, paste0("./lassosum2/prs/results/",race,"_sl_models.rds"))

res


