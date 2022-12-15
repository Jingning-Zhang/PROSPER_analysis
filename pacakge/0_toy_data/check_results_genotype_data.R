

# update geno snp id
library(bigreadr)
library(dplyr)
library(readr)
library(stringr)

for(eth in c("EUR","EAS","SAS","AFR","AMR")){

  bim <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",eth,"/tv_allchr.bim"))
  bim$RSID <- unlist(lapply(str_split(bim$V2,":"), FUN=function (x){x[1]}))
  tmp <- bim[,c("V2","RSID")]
  write_tsv(tmp, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/tmp/check_results_",eth,"_update_rsid"), col_names=F)

  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_sample_data/",eth,""))
  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 50 ",
  " --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",eth,"/tv_allchr ",
  " --keep /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",eth,"/validation.id.txt ",
  " --update-name /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/tmp/check_results_",eth,"_update_rsid ",
  " --make-bed ",
  " --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_sample_data/",eth,"/tuning_geno"))

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 50 ",
  " --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",eth,"/tv_allchr ",
  " --keep /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",eth,"/test.id.txt ",
  " --update-name /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/tmp/check_results_",eth,"_update_rsid ",
  " --make-bed ",
  " --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_sample_data/",eth,"/testing_geno"))

}

