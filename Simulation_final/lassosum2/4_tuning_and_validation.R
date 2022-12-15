rm(list = ls())

args <- commandArgs(T)
for(ii in 1:length(args)){ eval(parse(text=args[[ii]])) }
rep <- as.integer(rep)

library(bigreadr)
library(readr)

L=20
S=10

NCORES=12

system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
                " --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/", ethnic,"/tv_allchr",
                " --score /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/beta_rho_",rho,"_size_",size,".txt",
                " cols=+scoresums,-scoreavgs --score-col-nums 3-",L*S+2,
                " --out /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_",size))
SCORE_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_",size,".sscore"))
pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",ethnic, "/phenotypes_rho",rho,"_",setting,".phen"))

################################################
################################################

## tuning

tuning_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/test.id.txt"), stringsAsFactors = F)$V1
SCORE <- SCORE_all[match(tuning_id, SCORE_all$IID), -1:-4]

pheno <- pheno_all[match(tuning_id, pheno_all$V2), ]
pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]

R2 <- numeric(length=L*S)
for (i in 1:(L*S)){
  score <- scale(SCORE[,i], center = T, scale = T)
  if(sum(is.nan(score))>0){R2[i] <- 0; next}
  fit <- lm(pheno ~ score)
  R2[i] <- (coefficients(fit)['score'])^2/var(pheno)
}
best <- which.max(R2)

params2 <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/param_rho_",rho,"_size_",size,".txt"))
#params2[best,]

save(params2, best, file = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/best_tuning_rho_",rho,"_size_",size,".RData"))
saveRDS(R2, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/tuning_R2_rho_",rho,"_size_",size,".rds"))


################################################
################################################

## validation

validation_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/validation.id.txt"), stringsAsFactors = F)$V1
SCORE <- SCORE_all[match(validation_id, SCORE_all$IID), -1:-4]


pheno <- pheno_all[match(validation_id, pheno_all$V2), ]
pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep]

dat <- data.frame(y=as.numeric(pheno), score=SCORE[,best])
fit <- lm(y~., data = dat)
R2 <- summary(fit)$r.squared

saveRDS(R2, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,".rds"))

R2
