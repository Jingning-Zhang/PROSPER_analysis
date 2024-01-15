
library(bigreadr)

ethnic='SAS'

  print(ethnic)

  for(chr in 1:22){
    #print(chr)

    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block_UKBB/",ethnic,"/standard_data/chr",chr,"_snps.RData"))
    a <- unlist(snps_list)
    b <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/1000G/GRCh37/",ethnic,"/chr",chr,".bim"))
    b <- b[match(a, b$V2),]
    if(chr==1){
      res <- b
    }else{
      res <- rbind(res, b)
    }
    print(nrow(res))
    #if("" %in% b$V2){print(chr)}
  }
  fwrite2(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref/",ethnic,"/ref_bim.txt"), col.names = F, sep="\t", nThread=1)

# tar -czvf AFR_ukb.tar.gz AFR


res <- NULL
for(ethnic in c("AFR","AMR","EAS","EUR","SAS")){
  if(ethnic != "AMR"){
    dat <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref/",ethnic,"/ref_bim.txt"))
  }else{
    dat <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/1000g_ref/",ethnic,"/ref_bim.txt"))
  }

  if("" %in% dat$V2){print(ethnic)}

  res <- rbind(res,dat)
  res <- unique(res)
}
a = table(res$V2); max(a)
#res <- res[-nrow(res),]
fwrite2(res, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/ukbb_ref/ref_bim.txt"), col.names = F, sep="\t", nThread=1)
## this is same as the 1000g ref_bim.txt in /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/PROSPER/1000g_ref/ref_bim.txt
