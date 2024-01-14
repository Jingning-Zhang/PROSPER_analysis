
rm(list=ls())

library(bigsnpr)
library(bigreadr)
library(bigparallelr)

parameters = expand.grid(race = c('AFR','AMR','EUR'),
                         trait = c('height', 'bmi'))

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = as.character(parameters[as.numeric(temp1),1])
trait = as.character(parameters[as.numeric(temp1),2])

print(trait)
print(race)

NCORES <-  1

setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/AoU/",trait))

nsnps <- integer()
for (chr in 1:22){
  prs.file <- fread2(paste0("./ldpred2/coef/",race,"_ldpred2-chr",chr,".txt"))
  nsnps[chr] <- ncol(prs.file)
}

################################################
################################################

## tuning

for (chr in 1:22){

  print(paste0("tuning-",chr))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/",race,"/chr", chr,
                " --score ./ldpred2/coef/",race,"_ldpred2-chr",chr,".txt",
                " cols=+scoresums,-scoreavgs --score-col-nums 3-", nsnps[chr],
                " --out ./ldpred2/prs/tuning/",race,"_chr", chr,"_tuning_score"))
  score <- fread2(paste0("./ldpred2/prs/tuning/",race,"_chr", chr,"_tuning_score.sscore"))

  if(chr==1){
    SCORE_all <- score[, -1:-4,drop=F]
    SCORE_id <- score[,1:2,drop=F]
  }else{
    SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
  }
}
SCORE <- SCORE_all
save(SCORE, SCORE_id, file = paste0("./ldpred2/prs/tuning/",race,"_tuning_score.RData"))


################################################
################################################

## validation

dir.create("./ldpred2/prs/validation")

for (chr in 1:22){

  print(paste0("validation-",chr))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/",race,"/chr", chr,
                " --score ./ldpred2/coef/",race,"_ldpred2-chr",chr,".txt",
                " cols=+scoresums,-scoreavgs --score-col-nums 3-",nsnps[chr],
                " --out ./ldpred2/prs/validation/",race,"_chr", chr,"_validation_score"))

  score <- fread2(paste0("./ldpred2/prs/validation/",race,"_chr", chr,"_validation_score.sscore"))
  if(chr==1){
    SCORE_all <- score[, -1:-4,drop=F]
    SCORE_id <- score[,1:2,drop=F]
  }else{
    SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
  }
}
SCORE <- SCORE_all
save(SCORE, SCORE_id, file = paste0("./ldpred2/prs/validation/",race,"_validation_score.RData"))
