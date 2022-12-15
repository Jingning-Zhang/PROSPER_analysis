

# update geno snp id
library(bigreadr)
library(dplyr)
library(readr)
library(stringr)

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_5pop_lassosum2_param/lassosum2/")
for(eth in c("EUR","EAS","SAS","AFR","AMR")){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_5pop_lassosum2_param/lassosum2/",eth))
  load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_1_lassosum2/",eth,"/rep_1/best_tuning_rho_1_size_",ifelse(eth=="EUR",4,1),".RData"))
  a <- params2[best,c("delta","lambda")]
  colnames(a) <- c("delta0","lambda0")
  fwrite2(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/check_results_5pop_lassosum2_param/lassosum2/",eth,"/optimal_param.txt"), col.names = T, sep="\t", nThread=1)
}

