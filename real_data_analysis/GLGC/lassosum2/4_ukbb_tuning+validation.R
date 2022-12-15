 # /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

ethnic="AFR"
trait="TC"

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

path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2")

#############################################################
#############################################################

df <- fread2(paste0(path,"/Results/beta_",ethnic,"_",trait,".txt"), header = F, sep="\t", nThread=NCORES)
snps <- df$V1
sumdataA1 <- df$V2
sumdataA2 <- df$V3
beta_lassosum2 <- as.matrix(df[,-1:-3])

bim <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/",ethnic,"/allchr_bim"))

m <- which(!(snps %in% bim$V2))
bim <- bim[match(snps, bim$V2),]

if(length(m)>0) {
  snps <- snps[-m]
  sumdataA1 <- sumdataA1[-m]
  sumdataA2 <- sumdataA2[-m]
  bim <- bim[-m,]
  beta_lassosum2 <- beta_lassosum2[-m,]
}

qc <- allele.qc(a1=sumdataA1,a2=sumdataA2,ref1=bim$V5,ref2=bim$V6)

if(sum(qc$flip)>0){
  beta_lassosum2[qc$flip,] <- -1 * beta_lassosum2[qc$flip,]
}

df <- data.frame(SNP = snps, A1= bim$V5, beta_lassosum2, stringsAsFactors=F)
fwrite2(df, paste0(path,"/Results/ukbb_beta_",ethnic,"_",trait,".txt"), col.names = F, sep="\t", nThread=NCORES)

alltuning <- ncol(beta_lassosum2)

rm(list=c("df","beta_lassosum2","qc","bim","sumdataA1","sumdataA2","snps"))

#############################################################
#############################################################

#a <- read.table(paste0(path,"/Results/beta_",ethnic,"_",trait,".txt"), nrows=1)
#alltuning <- ncol(a)-2

dir.create(paste0(path,"/Results/test_all_prs/"))
dir.create(paste0(path,"/Results/final_results/"))

################################################
################################################

## tuning

for (chr in 1:22){

  print(paste0("tuning-",chr))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/tuning/",ethnic,"/chr", chr,
                " --score ", path,"/Results/ukbb_beta_",ethnic,"_",trait,".txt",
                " cols=+scoresums,-scoreavgs --score-col-nums 3-",alltuning+2,
                " --out ",path,"/Results/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_tuning_score"))

  score <- fread2(paste0(path,"/Results/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_tuning_score.sscore"))
  if(chr==1){
    SCORE_all <- score[, -1:-4,drop=F]
    SCORE_id <- score[,1:2,drop=F]
  }else{
    SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
  }
}
SCORE <- SCORE_all
save(SCORE, SCORE_id, file = paste0(path,"/Results/test_all_prs/",ethnic,"_",trait,"_tuning_score.RData"))


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

R2 <- numeric(length = alltuning)
for (i in 1:(alltuning)){
  model = lm(y~.-id,data = pheno); residual = model$residuals
  fit <- lm( residual ~ SCORE[,i] )
  R2[i] <- summary(fit)$r.square
}
best <- which.max(R2)
saveRDS(best, paste0(path,"/Results/final_results/",ethnic,"_",trait,"_best.rds"))


################################################
################################################

## validation

for (chr in 1:22){

  print(paste0("validation-",chr))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/genotype/validation/",ethnic,"/chr", chr,
                " --score ", path,"/Results/ukbb_beta_",ethnic,"_",trait,".txt",
                " cols=+scoresums,-scoreavgs --score-col-nums 3-",alltuning+2,
                " --out ",path,"/Results/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_validation_score"))

  score <- fread2(paste0(path,"/Results/test_all_prs/",ethnic,"_",trait,"_chr", chr,"_validation_score.sscore"))
  if(chr==1){
    SCORE_all <- score[, -1:-4,drop=F]
    SCORE_id <- score[,1:2,drop=F]
  }else{
    SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
  }
}
SCORE <- SCORE_all
save(SCORE, SCORE_id, file = paste0(path,"/Results/test_all_prs/",ethnic,"_",trait,"_validation_score.RData"))


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

score <- SCORE_validation[,best]

model = lm(y~.-id,data = pheno); residual = model$residuals
fit <- lm( residual~score )
R2 <- summary(fit)$r.square

tmp <- data.frame(best_R2=R2)
res <- tmp

saveRDS(res, paste0(path,"/Results/final_results/",ethnic,"_",trait,"_validation_R2_best.rds"))
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

saveRDS(res, paste0(path,"/Results/final_results/",ethnic,"_",trait,"_validation_R2_sl.rds"))

res_sl <- list(lasso=sl1, sl=sl2)
saveRDS(res_sl, paste0(path,"/Results/final_results/",ethnic,"_",trait,"_sl_models.rds"))

res

