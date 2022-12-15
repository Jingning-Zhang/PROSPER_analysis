

library(bigreadr)

path_to_sum = "/dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prscsx/sumdata/"

dir.create(path_to_sum)
snps <- NULL
for(eth in c("AFR","EUR","AMR","EAS","SAS")){
  sumdata <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/summdata/",eth,".txt"))
  sumdata$P <- 2*(1 - pnorm(abs(sumdata$beta/sumdata$beta_se)))
  sumdata <- sumdata[,c(1,3,4,5,8)]
  colnames(sumdata) <- c("SNP","A1","A2","BETA","P")
  fwrite2(sumdata, paste0(path_to_sum,"/",eth,".txt"), col.names = T, sep="\t", nThread=1)

  snps <- rbind(snps, fread2("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno.bim"))
  snps <- unique(snps)
}
fwrite2(snps, paste0(path_to_sum,"/alleth.bim"), col.names = F, sep="\t", nThread=1)

