#
#rm(list=ls())
#
#eid_geno <- bigreadr::fread2("/dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr1_v3.fam")$V1
#eid_geno_extract <- readLines("/dcl01/chatterj/data/ukbiobank/genetic/unrelated_european_eids/unrelated_european_ancestry.txt")[-1]
#eid_geno_extract <- eid_geno_extract[eid_geno_extract %in% eid_geno]
#eid_exclude <- readLines("/dcl01/chatterj/data/ukbiobank/phenotype/data/withdrawal_participants_UKBB_17712.csv")[-1]
#eid_geno_extract <- eid_geno_extract[!(eid_geno_extract %in% eid_exclude)]
#
#cov <- readRDS("/dcl01/chatterj/data/ukbiobank/phenotype/data/ukbiobank_covariates.rds")
#pheno <- readRDS("/dcl01/chatterj/data/ukbiobank/phenotype/data/ukb_11165_27504_April092019.rds")
#cov <- cov[cov$eid %in% eid_geno_extract, ]
#cov <- cov[!is.na(cov$eid),]
#
#eid <- cov$eid
#saveRDS(eid, "/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/ukbbeid.rds")
#
#set.seed(1)
#tmp <- sample(eid,5000)
#eid_tuning <- tmp[1:2500]
#eid_validation <- tmp[2501:5000]
#write.table(data.frame(eid_tuning,eid_tuning),
#            quote=F, row.names=F, col.names=F,
#            "/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/tuning.eid")
#write.table(data.frame(eid_validation,eid_validation),
#            quote=F, row.names=F, col.names=F,
#            "/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/validation.eid")
#
#
#cov_tuning <- cov[match(eid_tuning, cov$eid), c('eid','genetic_sex_f22001_0_0','age_when_attended_assessment_centre_f21003_0_0',paste0('genetic_principal_components_f22009_0_',1:10))]
#pheno_tuning <- pheno[match(eid_tuning, pheno$eid),c("eid","standing_height_f50_0_0")]
#save(cov_tuning, pheno_tuning, file = "/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/tuning.pheno")
#
#cov_validation <- cov[match(eid_validation, cov$eid), c('eid','genetic_sex_f22001_0_0','age_when_attended_assessment_centre_f21003_0_0',paste0('genetic_principal_components_f22009_0_',1:10))]
#pheno_validation <- pheno[match(eid_validation, pheno$eid),c("eid","standing_height_f50_0_0")]
#save(cov_validation, pheno_validation, file = "/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/validation.pheno")
#
##library(stringr)
##colnames(pheno)[str_detect(colnames(pheno),"height")]
#
#
#################################################################################################
#################################################################################################
#
#### genotype
#
#tmp <- paste0("#!/usr/bin/env bash
##$ -N tuning
##$ -l mem_free=1G,h_vmem=1G,h_fsize=100G
##$ -pe local 8
##$ -t 1-22
##$ -cwd
##$ -m e
##$ -M jzhan218@jhu.edu
#
#/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2",
#              " --threads 8",
#              " --extract /dcs04/nilanjan/data/jzhang2/MEPRS/mega-hm3-rsid.txt",
#                " --keep /dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/tuning.eid",
#                " --bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3",
#                " --rm-dup exclude-all",
#                " --make-bed",
#                " --out /dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/geno/tuning/chr$SGE_TASK_ID")
#
#writeLines(tmp, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/codes/tuning.sh"))
#
#
#
#tmp <- paste0("#!/usr/bin/env bash
##$ -N validation
##$ -l mem_free=1G,h_vmem=1G,h_fsize=100G
##$ -pe local 8
##$ -t 1-22
##$ -cwd
##$ -m e
##$ -M jzhan218@jhu.edu
#
#/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2",
#              " --threads 8",
#              " --extract /dcs04/nilanjan/data/jzhang2/MEPRS/mega-hm3-rsid.txt",
#                " --keep /dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/validation.eid",
#                " --bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3",
#                " --rm-dup exclude-all",
#                " --make-bed",
#                " --out /dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/geno/validation/chr$SGE_TASK_ID")
#
#writeLines(tmp, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning/codes/validation.sh"))
#

################################################################################################
################################################################################################
### SCORE

## EUR
load(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/EUR_height.RData"))

SCORE <- cbind(df_beta[,c(5,4)],beta_lassosum2)
write.table(SCORE,
            quote=F, row.names=F, col.names=F,
            "/dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum2.score")

#a <- bigreadr::fread2("/dcs04/nilanjan/data/23andme/Analysis/JZ/RunPRS/1_SingleEthnic/lassosum/EUR/height/prs.file")
#
#load("/dcs04/nilanjan/data/23andme/cleaned/snpinfo/snpinfo_mega.RData")
#a$V1 <- snpinfo_mega$assay.name[match(a$V1,snpinfo_mega$im.data.id)]
#sumdata <- bigreadr::fread2("/dcs04/nilanjan/data/23andme/cleaned/EUR/sumdat/height_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt")
#tmp <- sumdata$A1[match(a$V1,sumdata$rsid)]
#b <- cbind(a[,1],tmp,a[2:201])
#write.table(b,
#            quote=F, row.names=F, col.names=F,
#            "/dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum.score")


## AFR
load(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/AFR_height.RData"))

SCORE <- cbind(df_beta[,c(5,4)],beta_lassosum2)
write.table(SCORE,
            quote=F, row.names=F, col.names=F,
            "/dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/AFR_height_lassosum2.score")

#a <- bigreadr::fread2("/dcs04/nilanjan/data/23andme/Analysis/JZ/RunPRS/1_SingleEthnic/lassosum/AFR/height/prs.file")
#
#load("/dcs04/nilanjan/data/23andme/cleaned/snpinfo/snpinfo_mega.RData")
#a$V1 <- snpinfo_mega$assay.name[match(a$V1,snpinfo_mega$im.data.id)]
#sumdata <- bigreadr::fread2("/dcs04/nilanjan/data/23andme/cleaned/AFR/sumdat/height_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt")
#tmp <- sumdata$A1[match(a$V1,sumdata$rsid)]
#b <- cbind(a[,1],tmp,a[2:201])
#write.table(b,
#            quote=F, row.names=F, col.names=F,
#            "/dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/AFR_height_lassosum.score")


