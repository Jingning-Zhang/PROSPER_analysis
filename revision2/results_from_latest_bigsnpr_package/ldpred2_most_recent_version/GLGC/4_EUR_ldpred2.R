
rm(list=ls())

library(bigreadr)
library(readr)
library(dplyr)

parameters = expand.grid(race = c('AFR','AMR','EAS','EUR','SAS'),
                         trait = c('HDL','LDL','logTG','TC'))

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = as.character(parameters[as.numeric(temp1),1])
trait = as.character(parameters[as.numeric(temp1),2])

print(trait)
print(race)

NCORES <-  1

setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/GLGC/",trait))

dir.create(paste0("./ldpred2/prs/results/"))

################################################
################################################

## compute EUR prs

best = readRDS(paste0("./ldpred2/prs/tuning/EUR_best.rds"))

dir.create("./ldpred2/prs/validation")

for (chr in 1:22){

  print(paste0("validation-",chr))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/",race,"/chr", chr,
                " --score ./ldpred2/coef/EUR_ldpred2-chr",chr,".txt",
                " cols=+scoresums,-scoreavgs --score-col-nums ",best+2,
                " --out ./ldpred2/prs/validation/",race,"_chr", chr,"_validation_score_EUR"))

  score <- fread2(paste0("./ldpred2/prs/validation/",race,"_chr", chr,"_validation_score_EUR.sscore"))
  if(chr==1){
    SCORE_all <- score[, -1:-4,drop=F]
    SCORE_id <- score[,1:2,drop=F]
  }else{
    SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
  }
}
SCORE <- SCORE_all
save(SCORE, SCORE_id, file = paste0("./ldpred2/prs/validation/",race,"_validation_score_EUR.RData"))


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
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
pheno_validation <- pheno
SCORE_validation <- SCORE

score <- SCORE_validation

model = lm(y~.-id,data = pheno); residual = model$residuals
fit <- lm( residual~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(best_R2=R2)
res <- tmp

saveRDS(res, paste0("./ldpred2/prs/results/",race,"_validation_R2_EUR_best.rds"))
res

