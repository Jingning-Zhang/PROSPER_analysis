
rm(list=ls())
library(bigsnpr)
library(bigreadr)

races = c("AFR","AMR","EAS","EUR","SAS")

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = races[as.numeric(temp1)]
chr = as.integer(temp2)

ldr = 3/1000
NCORES = 1

path = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/",race)

setwd(path)
system(paste0("rm -rf ./ref_chr",chr,".bk"))
system(paste0("rm -rf ./ref_chr",chr,".rds"))

snp_readBed(paste0("./ref_chr",chr,".bed"))
obj.bigSNP <- snp_attach(paste0("./ref_chr",chr,".rds"))

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
map_ldref <- obj.bigSNP$map

af <- fread2(paste0("./freq_ref_chr",chr,".afreq"))
map_ldref$af <- af$ALT_FREQS[match(map_ldref$marker.ID, af$ID)]
map_ldref <- map_ldref[,-3]
colnames(map_ldref) <- c("chr","rsid","pos","a1","a0","a1_af")
saveRDS(map_ldref, paste0("./map_ldref_chr",chr,".rds"))

dir.create("tmp-data")
# To convert physical positions (in bp) to genetic positions (in cM), use
POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)

corr0 <- snp_cor(G, size = 3 / 1000, infos.pos = POS2, ncores = NCORES)
saveRDS(corr0, paste0("./LD_ref_chr",chr,".rds"))

