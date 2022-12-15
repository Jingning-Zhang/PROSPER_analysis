

# update geno snp id
library(bigreadr)
library(dplyr)
library(readr)
library(stringr)

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2_from_original_study/")
for(eth in c("EUR","AFR","AMR")){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2_from_original_study/",eth))
  best <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/final_results/",eth,"_height_best.rds"))
  params2 <- read_tsv(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/lassosum2/Results/para_",eth,"_height.txt"))
  a <- params2[best,c("delta","lambda")]
  colnames(a) <- c("delta0","lambda0")
  fwrite2(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_allofus_height/lassosum2_from_original_study/",eth,"/optimal_param.txt"), col.names = T, sep="\t", nThread=1)
}

