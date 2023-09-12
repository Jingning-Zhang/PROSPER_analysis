rm(list = ls())

args <- commandArgs(T)
for(ii in 1:length(args)){ eval(parse(text=args[[ii]])) }
rep <- as.integer(rep)

library(bigreadr)
library(readr)
library(SuperLearner)

#L=20
#S=10

#NCORES=12
#
#system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads ",NCORES,
#                " --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/", ethnic,"/tv_allchr",
#                " --score /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/beta_rho_",rho,"_size_",size,".txt",
#                " cols=+scoresums,-scoreavgs --score-col-nums 3-",L*S+2,
#                " --out /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_",size))
SCORE_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/score_rho_",rho,"_size_",size,".sscore"))
pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",ethnic, "/phenotypes_rho",rho,"_",setting,".phen"))


################################################
################################################

## super learning

## tuning

tuning_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/test.id.txt"), stringsAsFactors = F)$V1
SCORE <- SCORE_all[match(tuning_id, SCORE_all$IID), -1:-4]

pheno <- pheno_all[match(tuning_id, pheno_all$V2), ]
pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep,drop=F]

# Remove constant scores (marked by score_drop)
score_sd <- apply(SCORE,2,sd)
score_drop <- which(is.na(score_sd) | score_sd==0)
if(length(score_drop)>0){ SCORE <- SCORE[,-score_drop,drop=F] }


set.seed(20230814)
  sl = SuperLearner(Y = pheno,
                    X = data.frame(SCORE),
                    family = gaussian(),
                    SL.library = c("SL.glmnet", "SL.ridge", "SL.nnet","SL.lm") )


## validation

validation_id <- read.table(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/validation.id.txt"), stringsAsFactors = F)$V1
SCORE <- SCORE_all[match(validation_id, SCORE_all$IID), -1:-4]
pheno <- pheno_all[match(validation_id, pheno_all$V2), ]
pheno <- as.matrix(pheno[,-1:-2]); pheno <- pheno[,rep]

if(length(score_drop)>0){ SCORE <- SCORE[,-score_drop,drop=F] }

fit = lm(pheno ~ predict(sl, data.frame(SCORE), onlySL = TRUE)[[1]])
R2 <- summary(fit)$r.square

saveRDS(R2, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep,"/validation_R2_rho_",rho,"_size_",size,"_sl.rds"))

R2


