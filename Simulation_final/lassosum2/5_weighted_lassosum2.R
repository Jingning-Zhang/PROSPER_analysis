# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
rep <- as.integer(rep)

#setting=2
#rep=1
#ethnic="SAS"

library(bigreadr)
library(readr)
library(dplyr)


NCORES = 1

ethnic1_all <- setdiff(c("EAS","SAS","AMR","AFR"), ethnic)

for (rho in 1:3){
  for (size in 1:4){

    ################################################
    ## calculate PRS for non-target group
    for(ethnic1 in ethnic1_all ){
      load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/", ethnic1,"/rep_",rep,"/best_tuning_rho_",rho,"_size_",size,".RData"))

      system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                      " --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/", ethnic,"/tv_allchr",
                      " --score /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/", ethnic1,"/rep_",rep,"/beta_rho_",rho,"_size_",size,".txt",
                      " cols=+scoresums,-scoreavgs --score-col-nums ",best+2,
                      " --out /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_",size,"_", ethnic1,"_best"))
    }

    ################################################
    ## prepare data

    pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",ethnic, "/phenotypes_rho",rho,"_",setting,".phen"))
    SCORE_all_EUR <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_4_EUR_best.sscore"))
    SCORE_all_target <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_",size,".sscore"))
    SCORE_all_rest <- list()
    for(i in 1:length(ethnic1_all) ){
      ethnic1 <- ethnic1_all[i]
      SCORE_all_rest[[i]] <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_",size,"_", ethnic1,"_best.sscore"))
    }

    ################################################
    ## tuning

    tuning_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/test.id.txt"), stringsAsFactors = F)$V1
    pheno <- pheno_all[match(tuning_id, pheno_all$V2), ]
    pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]

    ## lassosum2 tuning score: target ethnic group
    SCORE <- SCORE_all_target[match(tuning_id, SCORE_all_target$IID), -1:-4]
    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/best_tuning_rho_",rho,"_size_",size,".RData"))
    score1 <- SCORE[,best]

    ## lassosum2 tuning score: EUR
    score2 <- SCORE_all_EUR[match(tuning_id, SCORE_all_EUR$IID), -1:-4]

    SCORE_final <- cbind(score1, score2)

    ## add all the other ancestry
    for(i in 1:length(ethnic1_all) ){
      ethnic1 <- ethnic1_all[i]
      tmp <- SCORE_all_rest[[i]]
      score_tmp <- tmp[match(tuning_id, tmp$IID), -1:-4]
      SCORE_final <- cbind(SCORE_final, score_tmp)
    }
    colnames(SCORE_final) <- c(ethnic, "EUR", ethnic1_all)


    ## weighted EUR+target
    tmp <- data.frame(pheno=as.numeric(pheno), SCORE_final[,c(ethnic, "EUR")])
    fit1 <- lm( pheno ~. , tmp )
    coeff1 <- coef(fit1)

    ## weighted all ancestry
    tmp <- data.frame(pheno=as.numeric(pheno), SCORE_final)
    fit2 <- lm( pheno ~ . , tmp )
    coeff2 <- coef(fit2)

    rm(list = c("pheno","tuning_id","SCORE","score1","score2","fit1","fit2","score_tmp","tmp","SCORE_final"))

    ################################################
    ## validation

    validation_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/validation.id.txt"), stringsAsFactors = F)$V1
    pheno <- pheno_all[match(validation_id, pheno_all$V2), ]
    pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]
    rm(list = c("pheno_all"))

    ## lassosum2 tuning score: target ethnic group
    SCORE <- SCORE_all_target[match(validation_id, SCORE_all_target$IID), -1:-4]
    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/best_tuning_rho_",rho,"_size_",size,".RData"))
    score1 <- SCORE[,best]

    ## lassosum2 tuning score: EUR
    score2 <- SCORE_all_EUR[match(validation_id, SCORE_all_EUR$IID), -1:-4]

    SCORE_final <- cbind(score1, score2)

    ## add all the other ancestry
    for(i in 1:length(ethnic1_all) ){
      ethnic1 <- ethnic1_all[i]
      tmp <- SCORE_all_rest[[i]]
      score_tmp <- tmp[match(validation_id, tmp$IID), -1:-4]
      SCORE_final <- cbind(SCORE_final, score_tmp)
    }
    colnames(SCORE_final) <- c(ethnic, "EUR", ethnic1_all)

    rm(list = c("SCORE_all_target", "SCORE_all_EUR", "SCORE_all_rest", "SCORE","best","validation_id","score_tmp","tmp"))


    ## weighted EUR+target
    tmp <- matrix(SCORE_final[,c(ethnic, "EUR")],ncol=2) %*% matrix(coeff1[c(ethnic, "EUR")],ncol = 1)
    fit <- lm(pheno~tmp)
    R2_2eth <- summary(fit)$r.squared

    ## weighted all ancestry
    tmp <- matrix(SCORE_final[,c(ethnic, "EUR", ethnic1_all)],ncol=5) %*% matrix(coeff2[c(ethnic, "EUR", ethnic1_all)],ncol = 1)
    fit <- lm(pheno~tmp)
    R2_5eth <- summary(fit)$r.squared

    R2 <- data.frame(R2_2eth=R2_2eth, R2_5eth=R2_5eth)

    saveRDS(R2, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,"_weighted_lassosum2.rds"))

    print(paste0(rho,"_",size))

    rm(list = c("pheno","score1","score2","coeff1","coeff2","fit","tmp","R2"))

  }
}

