


library(readr)
library(bigreadr)
library(dplyr)

mega <- readLines("/dcl01/chatterj/data/jin/MEGA/mega-hm3-rsid.txt")

for(race in c("EUR","AFR","AMR","EAS","SAS")){
  res <- tibble()
  for(chr in 1:22){
    print(chr)
    a <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/1000G/GRCh38/",race,"/chr",chr,"_freq.afreq"), col_types = cols())
    b <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/1000G/GRCh38/",race,"/chr",chr,".bim"), col_names=F, col_types = cols())
    a$POS <- b$X4[match(a$ID,b$X2)]
    res <- rbind(res, a[a$ID %in% mega, ])
  }
  write_tsv(res, paste0("/dcs04/nilanjan/data/jzhang2/1000G/GRCh38/",race,"/mega_freq.afreq") )
}








