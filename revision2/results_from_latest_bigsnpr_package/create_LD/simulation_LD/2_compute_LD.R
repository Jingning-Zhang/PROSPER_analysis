
rm(list=ls())

library(bigparallelr)
library(bigsnpr)
library(bigreadr)


races = c("AFR","AMR","EAS","EUR","SAS")

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = races[as.numeric(temp1)]
chr = as.integer(temp2)


# -------- GWAS summary statistics --------
# Summay statistics format

NCORES = 1

sumdata <- fread2(paste0("/dcl01/chatterj/data/hzhang1/multi_ethnic/LD_simulation_new/",race,"/pheno_summary_out_GA/summary_out_rho_1_size_1_rep_1_GA_1"), nThread = NCORES)
tmp <- stringr::str_split(sumdata$SNP,":")
sumdata$REF <- unlist(lapply(tmp, function (x){x[3]}))
sumdata$ALT <- unlist(lapply(tmp, function (x){x[4]}))
m <- sumdata$ALT != sumdata$A1
tmp <- sumdata$REF; tmp[m] <- sumdata$ALT[m]
sumdata$A0 <- tmp
sumdata$SE <- as.numeric(sumdata$BETA/sumdata$STAT)

sumstats <- sumdata[,c(1,2,3,12,4,7,13,9,6)]
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")

# -------- QC of summary statistics and reference data --------

# Match SNPs between summary statistics and reference data
ref <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/ref_allchr.bim"), nThread = NCORES)

snp <- intersect(ref$V2,sumstats$rsid)
ref <- ref[match(snp,ref$V2),]
sumstats <- sumstats[match(snp,sumstats$rsid),]
sumstats$pos <- ref$V4

flip <- ref$V5 != sumstats$a1
if(sum(flip)!=0){
  sumstats$beta[flip] <- -sumstats$beta[flip]
  sumstats$a1[flip] <- ref$V5[flip]
  sumstats$a0[flip] <- ref$V6[flip]
}
rm(list = c("ref","sumdata"))

# ------------------------ Reference Data preparation for LDpred2 ------------------------

tmp <- tryCatch({ rds <- snp_readBed(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/ref_allchr.bed")) }, error=function(e) e, warning=function(w) w)
if(is(tmp,"error")){ rds <- paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/ref_allchr.rds") }

obj.bigsnp <- snp_attach(rds)
# str(obj.bigsnp, max.level = 2, strict.width = "cut")

G   <- obj.bigsnp$genotypes
map <- dplyr::transmute(obj.bigsnp$map,
                        chr = chromosome, pos = physical.pos,
                        a0 = allele2, a1 = allele1, id = marker.ID)

info_snp <- snp_match(sumstats, map)
# Compute correlation between variants
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/LD_bigsnpr/",race))
workdir <- paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/LD_bigsnpr/",race)

#for (chr in 1:22) {

  print(chr)

  ## indices in 'sumstats'
  ind.chr <- which(info_snp$chr == chr)
  ## indices in 'G'
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

  POS2 <- snp_asGeneticPos(map$chr[ind.chr2], map$pos[ind.chr2], dir = workdir)
  tmp <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000, infos.pos = POS2, ncores = NCORES)
  saveRDS(tmp, paste0(workdir,"/corr_chr", chr, ".rds"))

  tmp <- map[ind.chr2,]
  saveRDS(tmp, paste0(workdir,"/map_chr", chr, ".rds"))

#}
#

