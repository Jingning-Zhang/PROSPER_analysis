
library(bigreadr)

#ethnic='AFR'
for(ethnic in c("EUR","SAS")){

  print(ethnic)

  for(chr in 1:22){
    print(chr)

    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block/",ethnic,"/standard_data/chr",chr,"_snps.RData"))
    a <- unlist(snps_list)
    b <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/1000G/GRCh37/",ethnic,"/chr",chr,".bim"))
    b <- b[match(a, b$V2),]
    if(chr==1){
      res <- b
    }else{
      res <- rbind(res, b)
    }
    print(nrow(res))
  }
  fwrite2(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data/",ethnic,"/ref_bim.txt"), col.names = F, sep="\t", nThread=1)

}

res <- NULL
for(ethnic in c("AFR","AMR","EAS","EUR","SAS")){
  dat <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data/",ethnic,"/ref_bim.txt"))
  res <- rbind(res,dat)
  res <- unique(res)
}
a = table(res$V2); max(a)
fwrite2(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/ref_data/ref_bim.txt"), col.names = F, sep="\t", nThread=1)

