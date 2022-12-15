# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

library(doMC)
library(foreach)


library(dplyr)
library(glmnet)

library(SuperLearner)
library(caret)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
L <- as.numeric(L); Lc <- as.numeric(Lc)
rep <- as.integer(rep)

NCORES <- 1
registerDoMC(NCORES)

#############################################################
#############################################################

M <- length(ethnic)
ethnic1 <- paste(ethnic,collapse = "_")

if(fixEUR=="T"){
  path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/4_multi-lassosum_by_chr_",para,"_fixEUR")
}else{
  path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para)
}

whichethnic <- which(ethnic == targetethnic)

alltuning <- integer(length(M))
for (mm in 1:M){
a <- system(paste0("awk -F'\t' '{print NF; exit}' ", path, "/Results/clean_results/Beta_rho_",rho,"_size_",size,"_",ethnic[mm],"_chr_1.txt"),
            intern = TRUE)
  alltuning[mm] <- as.integer(a)-2
}

for (mm in 1:M){

  for (chr in 1:22){

    print(paste0(mm,"-",chr))

    system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                  " --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/", targetethnic,"/tv_chr", chr,
                  " --score ", path,"/Results/clean_results/Beta_rho_",rho,"_size_",size,"_",ethnic[mm],"_chr_",chr,".txt",
                  " cols=+scoresums,-scoreavgs --score-col-nums 3", ifelse( (fixEUR == "T") & (ethnic[mm]=="EUR"), "", paste0("-",max(alltuning)+2)),
                  " --out ",path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_prs_rho_",rho,"_size_",size,"_chr", chr,"_tv_score"))

    score <- fread2(paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_prs_rho_",rho,"_size_",size,"_chr", chr,"_tv_score.sscore"))
    if(chr==1){
      SCORE_all <- score[, -1:-4,drop=F]
      SCORE_id <- score[,1:2,drop=F]
    }else{
      SCORE_all <- SCORE_all + score[, -1:-4,drop=F]
    }
  }
  SCORE <- SCORE_all
  save(SCORE, SCORE_id, file = paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_prs_rho_",rho,"_size_",size,"_tv_score.RData"))

}

for (mm in 1:M){
  load(paste0(path,"/Results/test_all_prs/",targetethnic,"/",ethnic[mm],"_prs_rho_",rho,"_size_",size,"_tv_score.RData"))
  if(mm==1){
    SCORE_all <- SCORE
    SCORE_id_all <- SCORE_id
  }else{
    SCORE_all <- cbind(SCORE_all, SCORE[match(SCORE_id_all$IID, SCORE_id$IID),])
  }
}
SCORE_id <- SCORE_id_all
save(SCORE_all, SCORE_id,
     file = paste0(path,"/Results/test_all_prs/",targetethnic,"/all_prs_rho_",rho,"_size_",size,"_tv_score.RData"))


################################################
pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",targetethnic, "/phenotypes_rho",rho,"_",setting,".phen"))


################################################
################################################

## best tuning prs

nprs <- ifelse(fixEUR == "T", (M-1)*max(alltuning)+1, M*max(alltuning))
prseth  <- character()
for (mm in 1:M){
  if( (fixEUR == "T") & (ethnic[mm]=="EUR") ){
    prseth <- c(prseth, "EUR")
  }else{
    prseth <- c(prseth, rep(ethnic[mm], max(alltuning)))
  }
}

################################################
################################################

## tuning

id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",targetethnic,"/test.id.txt"), stringsAsFactors = F)$V2
pheno <- pheno_all[match(id, pheno_all$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,,drop=F] }

R2 <- numeric(length = nprs)
for (i in 1:(nprs)){
  score <- scale(SCORE[,i], center = T, scale = T)
  if(sum(is.nan(score))>0){R2[i] <- 0; next}
  fit <- lm(pheno ~ score)
  R2[i] <- (coefficients(fit)['score'])^2/var(pheno)
}
best <- which.max(R2)
best_target <- which(prseth==targetethnic)[which.max(R2[prseth==targetethnic])]

## validation

id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",targetethnic,"/validation.id.txt"), stringsAsFactors = F)$V2
pheno <- pheno_all[match(id, pheno_all$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,,drop=F] }

dat <- data.frame(y=as.numeric(pheno), score=SCORE[,best_target])
fit <- lm(y~., data = dat)
R2 <- summary(fit)$r.squared

tmp <- data.frame(best_R2=R2)
res <- tmp

dat <- data.frame(y=as.numeric(pheno), score=SCORE[,best])
fit <- lm(y~., data = dat)
R2 <- summary(fit)$r.squared

tmp <- data.frame(best_R2_across_ancestry=R2)
res <- cbind(res,tmp)


################################################
################################################

## combined prs

################################################
################################################

## tuning

id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",targetethnic,"/test.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(id, pheno_all$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,,drop=F] }
pheno_tuning <- pheno
SCORE_tuning = SCORE


## validation
id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",targetethnic,"/validation.id.txt"), stringsAsFactors = F)$V1
pheno <- pheno_all[match(id, pheno_all$V2), ]; pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
SCORE <- SCORE_all[match(id, SCORE_id$IID),];  SCORE <- as.matrix(SCORE)
m <- unique(which(is.na(SCORE), arr.ind=TRUE)[,"row"]); if(length(m)>0){ SCORE <- SCORE[-m,]; pheno <- pheno[-m,,drop=F] }
pheno_validation <- pheno
SCORE_validation = SCORE
rm(list = c("pheno","SCORE"))


set.seed(20220502)
  sl = SuperLearner(Y = as.numeric(pheno_tuning),
                    X = data.frame(SCORE_tuning[,prseth==targetethnic,drop=F]),
                    family = gaussian(),
                    SL.library = c("SL.glmnet") )
  score = predict(sl, data.frame(SCORE_validation[,prseth==targetethnic,drop=F]), onlySL = TRUE)[[1]]
  fit <- lm(pheno_validation ~ score)
  R2 <- summary(fit)$r.square

tmp <- data.frame(lasso_combined_R2=R2)
res <- cbind(res,tmp)



set.seed(20220502)
  sl = SuperLearner(Y = as.numeric(pheno_tuning),
                    X = data.frame(SCORE_tuning),
                    family = gaussian(),
                    SL.library = c("SL.glmnet") )
  score = predict(sl, data.frame(SCORE_validation), onlySL = TRUE)[[1]]
  fit <- lm(as.numeric(pheno_validation) ~ score)
  R2 <- summary(fit)$r.square
tmp <- data.frame(lasso_combined_R2_cross_ancestry=R2)
res <- cbind(res,tmp)

saveRDS(res, paste0(path,"/Results/test_all_prs/",targetethnic,"/rho_",rho,"_size_",size,"_validation_R2_no_sl.rds"))
apply(res,MARGIN = 2, mean)


set.seed(20220502)
  sl = SuperLearner(Y = as.numeric(pheno_tuning),
                    X = data.frame(SCORE_tuning[,prseth==targetethnic,drop=F]),
                    family = gaussian(),
                    SL.library = c("SL.glmnet", "SL.ridge", "SL.nnet","SL.lm") )
  score = predict(sl, data.frame(SCORE_validation[,prseth==targetethnic,drop=F]), onlySL = TRUE)[[1]]
  fit <- lm(as.numeric(pheno_validation) ~ score)
  R2 <- summary(fit)$r.square
tmp <- data.frame(superlearning_combined_R2=R2)
res <- cbind(res,tmp)


set.seed(20220502)
  sl = SuperLearner(Y = as.numeric(pheno_tuning),
                    X = data.frame(SCORE_tuning),
                    family = gaussian(),
                    SL.library = c("SL.glmnet", "SL.ridge", "SL.nnet","SL.lm") )
  score = predict(sl, data.frame(SCORE_validation), onlySL = TRUE)[[1]]
  fit <- lm(as.numeric(pheno_validation) ~ score)
  R2 <- summary(fit)$r.square
tmp <- data.frame(superlearning_combined_R2_cross_ancestry=R2)
res <- cbind(res,tmp)

saveRDS(res, paste0(path,"/Results/test_all_prs/",targetethnic,"/rho_",rho,"_size_",size,"_validation_R2.rds"))

res


