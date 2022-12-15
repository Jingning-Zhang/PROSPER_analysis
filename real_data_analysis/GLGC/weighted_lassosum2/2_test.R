# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

ethnic="AFR"
trait="HDL"

rm(list=ls())

library(readr)
library(bigreadr)
library(dplyr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

NCORES <- 1

path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/weighted_lassosum2")

#############################################################
#############################################################

a <- read.table(paste0(path,"/ukbb_beta_",trait,".txt"), nrows=1)
alltuning <- ncol(a)-2



################################################
################################################

## tuning

for (chr in 1:22){

  print(paste0("tuning-",chr))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/",ethnic,"/chr", chr,
                " --score ", path,"/ukbb_beta_",trait,".txt",
                " header-read cols=+scoresums,-scoreavgs --score-col-nums 3-",alltuning+2,
                " --out ",path,"/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_tuning_score"))

  score <- fread2(paste0(path,"/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_tuning_score.sscore"))
  if(chr==1){
    SCORE_all <- score[, -1:-4,drop=F]
    SCORE_id <- score[,1:2,drop=F]
  }else{
    SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
  }
}
SCORE <- SCORE_all
save(SCORE, SCORE_id, file = paste0(path,"/test_all_prs/",ethnic,"_",trait,"_tuning_score.RData"))


pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",ethnic,"_tuning.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",ethnic,"_tuning.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

iid <- intersect(pheno$id, SCORE_id$IID)
pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
pheno_tuning <- pheno
SCORE_tuning <- SCORE


model = lm(y~.-id,data = pheno_tuning); residual = model$residuals
fit_2 <- lm( residual~., data.frame(residual,SCORE_tuning[,colnames(SCORE_tuning) %in% paste0(c("EUR",ethnic),"_SUM")]) )
fit_all <- lm( residual~., data.frame(residual,SCORE_tuning) )


################################################
################################################

## validation

for (chr in 1:22){

  print(paste0("validation-",chr))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/",ethnic,"/chr", chr,
                " --score ", path,"/ukbb_beta_",trait,".txt",
                " header-read cols=+scoresums,-scoreavgs --score-col-nums 3-",alltuning+2,
                " --out ",path,"/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_validation_score"))

  score <- fread2(paste0(path,"/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_validation_score.sscore"))
  if(chr==1){
    SCORE_all <- score[, -1:-4,drop=F]
    SCORE_id <- score[,1:2,drop=F]
  }else{
    SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
  }
}
SCORE <- SCORE_all
save(SCORE, SCORE_id, file = paste0(path,"/test_all_prs/",ethnic,"_",trait,"_validation_score.RData"))


pheno = read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/",trait,"/tuning+validation/",ethnic,"_validation.txt"))
pheno = pheno[,1:2]
covar <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/",ethnic,"_validation.txt"))
pheno <- left_join(pheno, covar)
colnames(pheno) = c('id','y','sex','age',paste0('pc',1:10))
pheno = pheno[complete.cases(pheno$y),]

iid <- intersect(pheno$id, SCORE_id$IID)
pheno <- pheno[match(iid, pheno$id), ]
SCORE <- as.matrix(SCORE)[match(iid, SCORE_id$IID), ]
pheno_validation <- pheno
SCORE_validation <- SCORE


model = lm(y~.-id,data = pheno_validation); residual = model$residuals

score = predict(fit_2, data.frame(SCORE_validation[,colnames(SCORE_validation) %in% paste0(c("EUR",ethnic),"_SUM")]))
fit <- lm( residual~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(weighted_2=R2)
res <- tmp

score = predict(fit_all, data.frame(SCORE_validation))
fit <- lm( residual~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(weighted_all=R2)
res <- cbind(res,tmp)

saveRDS(res, paste0(path,"/final_results/",ethnic,"_",trait,"_validation_R2_weighted.rds"))

res

