 rm(list = ls())

#ethnic="AFR"
#trait="HDL"

args <- commandArgs(T)
for(ii in 1:length(args)){ eval(parse(text=args[[ii]])) }

library(bigparallelr)
library(bigsnpr)
library(bigreadr)
library(readr)

# ----------------------- Functions -----------------------
# Allele match function
allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)

  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;

  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  snp[["amb_str"]] = (a1 == flip2 & a2 == flip1) | (a1 == flip1 & a2 == flip2)

  return(snp)
}

# -------- GWAS summary statistics --------
# Summay statistics format
# Colnames: chromosome <- "chr", base pair position <- "pos", snp rsid <- "rsid", effect allele <- "a1", other allele <- "a0",
#           beta <- "beta", standard error <- "beta_se", case sample size <- "n_case", control sample size <- "n_control", P value <- "p"

#NCORES <-  nb_cores()
NCORES = 1

sumdata <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",ethnic,"/",trait,".txt"))
#if(trait %in% c("any_cvd", "depression", "iqb.sing_back_musical_note", "migraine_diagnosis", "morning_person")){
#  sumdata$n_eff <- 4 / (1 / sumdata$N_control + 1 / sumdata$N_case)
#
#}else{
#  sumdata$n_eff <- sumdata$N
#}

sumstats <- sumdata[,c(2,1,3,10,9,5,6,7,4)]
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")

# -------- QC of summary statistics and reference data --------

# Match SNPs between summary statistics and reference data
ref <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/",ethnic,"/",trait,"/allchr.bim"))

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

tmp <- tryCatch({ rds <- snp_readBed(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/",ethnic,"/",trait,"/allchr.bed")) }, error=function(e) e, warning=function(w) w)
if(is(tmp,"error")){ rds <- paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/",ethnic,"/",trait,"/allchr.rds") }

obj.bigsnp <- snp_attach(rds)
# str(obj.bigsnp, max.level = 2, strict.width = "cut")

G   <- obj.bigsnp$genotypes
map <- dplyr::transmute(obj.bigsnp$map,
                        chr = chromosome, pos = physical.pos,
                        a0 = allele2, a1 = allele1)

CHR <- obj.bigsnp$map$chromosome
POS <- obj.bigsnp$map$physical.pos

info_snp <- snp_match(sumstats, map, strand_flip = FALSE)
# Compute correlation between variants
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/intermediate"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/intermediate/",ethnic))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/intermediate/",ethnic,"/",trait))
workdir <- paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/intermediate/",ethnic,"/",trait)
df_beta <- info_snp

#for (chr in 1:22) {
#
#  print(chr)
#  corr0 <- runonce::save_run({
#    ## indices in 'sumstats'
#    ind.chr <- which(df_beta$chr == chr)
#    ## indices in 'G'
#    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
#
#    POS2 <- snp_asGeneticPos(map$chr[ind.chr2], map$pos[ind.chr2], dir = workdir)
#    snp_cor(G, ind.col = ind.chr2, size = 3 / 1000, infos.pos = POS2,
#            ncores = NCORES)
#
#  }, file =  paste0(workdir,"/corr_chr", chr, ".rds"))
#
#}


system(paste0("rm ",workdir,"/corr_chr1_tmp.sbk"))

for (chr in 1:22) {

  print(chr)

  corr0 <- readRDS(paste0(workdir,"/corr_chr", chr, ".rds"))

  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, paste0(workdir,"/corr_chr", chr, "_tmp"), compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
file.size(corr$sbk) / 1024^3

print(paste0('Complete data preparation'))

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results"))
save(df_beta, corr, file = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/",ethnic,"_",trait,"_corr.RData"))
saveRDS(df_beta, file = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/",ethnic,"_",trait,"_df_beta.rds"))

