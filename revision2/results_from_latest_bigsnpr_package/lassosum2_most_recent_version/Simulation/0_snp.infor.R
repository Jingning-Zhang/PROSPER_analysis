
rm(list=ls())

load("/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/snp.infor.rdata")

pos <- integer()
for(race  in  c('AFR','AMR','EAS','EUR','SAS')){
  map_ldref <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/",race,"/map_ldref.rds"))
  pos <- c(pos,map_ldref$pos)
}

snp.infor1 <- snp.infor[snp.infor$position %in% pos,]
saveRDS(snp.infor1, "/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/snp.infor.rds")



