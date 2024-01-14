
rm(list=ls())
library(bigsnpr)
library(bigreadr)

races = c("AFR","AMR","EAS","EUR","SAS")

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = races[as.numeric(temp1)]


#for(race in races){
  print(race)

  workdir <- paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/LD_bigsnpr/",race)

  for (chr in 1:22) {

    print(chr)

    corr0 <- readRDS(paste0(workdir,"/corr_chr", chr, ".rds"))
    map0 <- readRDS(paste0(workdir,"/map_chr", chr, ".rds"))

    if (chr == 1) {
      system(paste0("rm ",workdir,"/corr_chr", chr, "_tmp.sbk"))
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, paste0(workdir,"/corr_chr", chr, "_tmp"), compact = TRUE)
      map <- map0
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
      map <- rbind(map, map0)
    }
  }
  file.size(corr$sbk) / 1024^3

  saveRDS(corr, paste0(workdir,"/corr_sfbm.rds"))
  saveRDS(map, paste0(workdir,"/map.rds"))
  saveRDS(ld, paste0(workdir,"/ld.rds"))


#}

