rm(list = ls())

args <- commandArgs(T)
for(ii in 1:length(args)){ eval(parse(text=args[[ii]])) }

library(bigparallelr)
library(bigsnpr)
library(bigreadr)
library(readr)

# -------- GWAS summary statistics --------
# Summay statistics format

NCORES = 2

sumdata <- fread2(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/",ethnic,"/pheno_summary_out_GA/summary_out_rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting),
                  nThread = NCORES)
tmp <- stringr::str_split(sumdata$SNP,":")
sumdata$REF <- unlist(lapply(tmp, function (x){x[3]}))
sumdata$ALT <- unlist(lapply(tmp, function (x){x[4]}))
m <- sumdata$ALT != sumdata$A1
tmp <- sumdata$REF; tmp[m] <- sumdata$ALT[m]
sumdata$A0 <- tmp
sumdata$SE <- as.numeric(sumdata$BETA/sumdata$STAT)

sumstats <- sumdata[,c(1,2,3,12,4,7,13,9,6)]
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
dir.create( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic) )
write_tsv(sumstats, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic,"/rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting,".txt"))


# -------- QC of summary statistics and reference data --------

# Match SNPs between summary statistics and reference data
ref <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/ref_allchr.bim"))

snp <- intersect(ref$V2,sumstats$rsid)
ref <- ref[match(snp,ref$V2),]
sumstats <- sumstats[match(snp,sumstats$rsid),]
sumstats$pos <- ref$V4

flip <- ref$V5 != sumstats$a1
if(sum(flip)!=0){
  sumstats$beta[flip] <- -sumstats$beta[flip]
  sumstats$a1[flip] <- ref$V5[flip]
  sumstats$a0[flip] <- ref$V6[flip]
}
rm(list = c("ref","sumdata"))

map <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr/",ethnic,"/map.rds"))
df_beta <- snp_match(sumstats, map, strand_flip = FALSE)
saveRDS(df_beta, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic,"/rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting,"_matched_to_ref.rds"))

betamax <- max( abs( sumstats$beta / sqrt(sumstats$n_eff * sumstats$beta_se^2 + sumstats$beta^2) ) )
for (chr in 1:22){
  map <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/LD_bigsnpr/",ethnic,"/map_chr", chr, ".rds"))
  df_beta <- snp_match(sumstats, map, strand_flip = FALSE)
  save(df_beta, betamax, file=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic,"/rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting,"_matched_to_ref_chr",chr,".RData"))
}

