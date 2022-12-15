
library(bigreadr)

pheno_all <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/AMR/phenotypes_rho1_1.phen"))[,1:3]
fwrite2(pheno_all, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AMR/pheno.fam"), col.names = F, sep=" ", nThread=1)





