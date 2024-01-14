
rm(list=ls())

library(bigreadr)
library(readr)
library(dplyr)

parameters = expand.grid(race = c('AFR','AMR'),
                         trait = c('height', 'bmi'))

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = as.character(parameters[as.numeric(temp1),1])
trait = as.character(parameters[as.numeric(temp1),2])

print(race)
print(trait)

NCORES <-  1

setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait))

dir.create(paste0("./lassosum2/prs/results/"))

################################################
################################################

## compute EUR prs

best = readRDS(paste0("./lassosum2/prs/tuning/EUR_best.rds"))

load(paste0("./lassosum2/prs/validation/",race,"-EUR_validation_score.RData"))


################################################
################################################

## testing performance

pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",race,"_validation.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",race,"_validation.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

iid <- intersect(pheno$id, SCORE_id$IID)
pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), best, drop=F]
pheno_validation <- pheno
SCORE_validation <- SCORE

score <- SCORE_validation

model = lm(y~.-id,data = pheno); residual = model$residuals
fit <- lm( residual~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(EUR_R2=R2)
res <- tmp

saveRDS(res, paste0("./lassosum2/prs/results/",race,"_validation_R2_EUR_best.rds"))
res

