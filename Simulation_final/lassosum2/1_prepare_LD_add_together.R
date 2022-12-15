
rm(list = ls())

args <- commandArgs(T)
for(ii in 1:length(args)){ eval(parse(text=args[[ii]])) }
library(bigparallelr)
library(bigsnpr)

ethnic = 'AFR'
# -------- GWAS summary statistics --------
# Summay statistics format

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr/",ethnic))
workdir <- paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr/",ethnic)

for (chr in 1:22) {

  print(chr)

  corr0 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/LD_bigsnpr/",ethnic,"/corr_chr", chr, ".rds"))
  map0 <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/LD_bigsnpr/",ethnic,"/map_chr", chr, ".rds"))

  if (chr == 1) {
    corr <- as_SFBM(corr0, paste0(workdir,"/corr_chr", chr, "_tmp"), compact = TRUE)
    map <- map0
  } else {
    corr$add_columns(corr0, nrow(corr))
    map <- rbind(map, map0)
  }
}
file.size(corr$sbk) / 1024^3

saveRDS(corr, paste0(workdir,"/corr_sfbm.rds"))
saveRDS(map, paste0(workdir,"/map.rds"))



