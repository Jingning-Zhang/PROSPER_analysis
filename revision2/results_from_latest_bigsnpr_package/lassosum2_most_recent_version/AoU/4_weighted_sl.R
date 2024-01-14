# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())


library(bigreadr)
library(readr)
library(dplyr)
library(SuperLearner)

parameters = expand.grid(race = c('AFR','AMR'),
                         trait = c('height', 'bmi'))

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = as.character(parameters[as.numeric(temp1),1])
trait = as.character(parameters[as.numeric(temp1),2])

print(race)
print(trait)


setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait))

dir.create(paste0("./lassosum2/prs/results/"))

NCORES = 1


################################################
## prepare data
################################################

## score

race1_all <- c('AFR','AMR','EUR')

## best tuning parameters
best <- NULL
for(i in 1:length(race1_all) ){
  race1 <- race1_all[i]
  best <- c(best, readRDS(paste0("./lassosum2/prs/tuning/",race1,"_best.rds")))
}
names(best) <- race1_all

################################################
## tuning data

pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",race,"_tuning.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",race,"_tuning.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

for(i in 1:length(race1_all) ){
  race1 <- race1_all[i]
  load(paste0("./lassosum2/prs/tuning/",race,"-",race1,"_tuning_score.RData"))
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

iid <- intersect(pheno$id, SCORE_id$IID)

pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
SCORE_best <- as.matrix(SCORE_best)[match(iid, SCORE_id$IID), ]

model = lm(y~.-id,data = pheno)
pheno_tuning <- model$residuals
SCORE_tuning <- SCORE
SCORE_tuning_best <- SCORE_best

################################################
## validation data

pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",race,"_validation.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",race,"_validation.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

for(i in 1:length(race1_all) ){
  race1 <- race1_all[i]
  load(paste0("./lassosum2/prs/validation/",race,"-",race1,"_validation_score.RData"))

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

iid <- intersect(pheno$id, SCORE_id$IID)

pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
SCORE_best <- as.matrix(SCORE_best)[match(iid, SCORE_id$IID), ]

model = lm(y~.-id,data = pheno)
pheno_validation <- model$residuals
SCORE_validation <- SCORE
SCORE_validation_best <- SCORE_best

################################################
## test


# super learning

tmp <- data.frame(SCORE_tuning); m <- apply(tmp,2,sd)==0; tmp <- tmp[,!m, drop=F]
sl = SuperLearner(Y = pheno_tuning,
                  X = tmp,
                  family = gaussian(),
                  SL.library = c("SL.glmnet", "SL.ridge","SL.lm") )

score <- predict(sl, data.frame(SCORE_validation)[,!m, drop=F], onlySL = TRUE)[[1]]

fit = lm(pheno_validation~score)
R2 <- summary(fit)$r.square
tmp <- data.frame(superlearning_combined_R2=R2)
res <- tmp


# weighted

tmp <- data.frame(pheno_tuning,SCORE_tuning_best); colnames(tmp) <- c("y",race1_all)
wt = lm(y~., tmp)
tmp <- data.frame(SCORE_validation_best); colnames(tmp) <- race1_all
score <- predict(wt, tmp)

fit = lm(pheno_validation~score)
R2 <- summary(fit)$r.square
tmp <- data.frame(weighted_combined_R2=R2)
res <- cbind(res,tmp)

# best

fit = lm(pheno_validation~SCORE_validation_best[,which(race1_all==race)])
R2 <- summary(fit)$r.square
tmp <- data.frame(best_R2=R2)
res <- cbind(res,tmp)

saveRDS(res, paste0("./lassosum2/prs/results/",race,"_validation_R2_sl.rds"))

res_model <- list(sl=sl, wt = wt)
saveRDS(res_model, paste0("./lassosum2/prs/results/",race,"_sl_models.rds"))

res

