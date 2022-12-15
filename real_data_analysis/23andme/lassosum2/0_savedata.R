rm(list = ls())

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

NCORES = 8

ethnic='AFR'


#for (trait in c("depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){ # EUR

for (trait in c("any_cvd","depression","heart_metabolic_disease_burden")){ # AFR

  print(trait)

  # -------- GWAS summary statistics --------

  sumdata <- bigreadr::fread2(paste0("/dcs04/nilanjan/data/23andme/cleaned/",ethnic,"/sumdat/",trait,"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"))
  if(trait %in% c("any_cvd", "depression", "iqb.sing_back_musical_note", "migraine_diagnosis", "morning_person")){
    sumdata$n_eff <- 4 / (1 / sumdata$N_control + 1 / sumdata$N_case)
  }else{
    sumdata$n_eff <- sumdata$N_control
  }

  sumstats <- sumdata[,c(2,1,3,5,4,10,11,9,12)]
  names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")

  # -------- QC of summary statistics and reference data --------

  # Match SNPs between summary statistics and reference data
  ref <- fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/ref_geno/",ethnic,"/",trait,"/allchr.bim"))

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
  rm(list = c("ref","sumdata","snp"))

  # ------------------------ Reference Data preparation for LDpred2 ------------------------

  rds <- paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/ref_geno/",ethnic,"/",trait,"/allchr.rds")
  obj.bigsnp <- snp_attach(rds)
  # str(obj.bigsnp, max.level = 2, strict.width = "cut")

  G   <- obj.bigsnp$genotypes
  map <- dplyr::transmute(obj.bigsnp$map,
                          chr = chromosome, pos = physical.pos,
                          a0 = allele2, a1 = allele1)

  info_snp <- snp_match(sumstats, map, strand_flip = FALSE)
  # Compute correlation between variants
  dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/",ethnic))
  dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/",ethnic,"/",trait))
  workdir <- paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/",ethnic,"/",trait)
  df_beta <- info_snp

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

  save(df_beta, corr, file = paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/",ethnic,"_",trait,"_corr.RData"))

  rm(list = c("df_beta","corr","corr0","ld","info_snp","G","map","obj.bigsnp","workdir"))

}
