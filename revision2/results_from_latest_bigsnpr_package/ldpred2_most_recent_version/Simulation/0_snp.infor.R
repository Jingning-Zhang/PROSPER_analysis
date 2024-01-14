
rm(list=ls())

load("/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/snp.infor.rdata")

pos <- integer()
for(race  in  c('AFR','AMR','EAS','EUR','SAS')){
  map_ldref <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/",race,"/map_ldref.rds"))
  pos <- c(pos,map_ldref$pos)
}

snp.infor1 <- snp.infor[snp.infor$position %in% pos,]
saveRDS(snp.infor1, "/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/snp.infor.rds")


sumdata <- fread2(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/",race,"/pheno_summary_out_GA/summary_out_rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting),
                  nThread = NCORES)

map_ldref <- readRDS( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/updated_bigsnpr/map_ldref.rds"))
sumdata <- sumdata[sumdata$SNP %in% map_ldref$rsid,]

sumdata$REF <- map_ldref$a0[match(sumdata$SNP, map_ldref$rsid)]
sumdata$ALT <- map_ldref$a1[match(sumdata$SNP, map_ldref$rsid)]
m <- sumdata$ALT != sumdata$A1
tmp <- sumdata$REF; tmp[m] <- sumdata$ALT[m]
sumdata$A0 <- tmp
sumdata$SE <- as.numeric(sumdata$BETA/sumdata$STAT)

sumstats <- sumdata[,c(1,2,3,12,4,7,13,9,6)]
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
dir.create( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic) )
write_tsv(sumstats, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic,"/rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting,".txt"))


