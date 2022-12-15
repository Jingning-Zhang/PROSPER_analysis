# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

#ethnic="AFR"
#trait="HDL"

rm(list=ls())

allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)

  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;

  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  snp[["amb_str"]] = (a1 == flip2 & a2 == flip1) | (a1 == flip1 & a2 == flip2)

  return(snp)
}

library(readr)
library(bigreadr)


library(dplyr)
library(glmnet)

library(SuperLearner)
library(caret)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

NCORES <- 1

path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2")

#############################################################
#############################################################

a <- read.table(paste0(path,"/Results/ukbb_beta_",ethnic,"_",trait,".txt"), nrows=1)
alltuning <- ncol(a)-2

################################################
################################################

## best EUR

best_EUR <- readRDS(paste0(path,"/Results/final_results/EUR_",trait,"_best.rds"))

################################################
################################################

## validation

for (chr in 1:22){

  print(paste0("validation-",chr))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/",ethnic,"/chr", chr,
                " --score ", path,"/Results/ukbb_beta_EUR_",trait,".txt",
                " cols=+scoresums,-scoreavgs --score-col-nums ",best_EUR+2,
                " --out ",path,"/Results/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_validation_score_best_EUR_lassosum2"))

  score <- fread2(paste0(path,"/Results/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_validation_score_best_EUR_lassosum2.sscore"))
  if(chr==1){
    SCORE_all <- score[, -1:-4,drop=F]
    SCORE_id <- score[,1:2,drop=F]
  }else{
    SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
  }
}
SCORE <- SCORE_all
save(SCORE, SCORE_id, file = paste0(path,"/Results/test_all_prs/",ethnic,"_",trait,"_validation_score_best_EUR_lassosum2.RData"))


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

model = lm(y~.-id,data = pheno); residual = model$residuals
fit <- lm( residual~SCORE )
R2 <- summary(fit)$r.square

tmp <- data.frame(best_EUR_R2=R2)
res <- tmp

saveRDS(res, paste0(path,"/Results/final_results/",ethnic,"_",trait,"_validation_R2_best_EUR_lassosum2.rds"))
res


