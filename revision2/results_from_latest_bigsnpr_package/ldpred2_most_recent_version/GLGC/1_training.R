
rm(list=ls())

library(bigsnpr)
library(bigreadr)
library(bigparallelr)

races = c('AFR','AMR','EAS','EUR','SAS')
traits = c('HDL','LDL','logTG','TC')

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = races[as.numeric(temp1)]
trait = traits[as.numeric(temp2)]

print(trait)
print(race)

NCORES <- 3

map_ldref <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/",race,"/map_ldref.rds"))
sumraw = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",race,"/",trait,".txt"))
sumstats = sumraw[,c('CHR','rsID','POS_b37','A2','A1','BETA','SE','P','N','A1_FREQ_1000G')]
names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff","a1_sumdata_af")

info_snp <- snp_match(sumstats, map_ldref, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
info_snp <- tidyr:: drop_na(tibble::as_tibble(info_snp))
sd_ldref <- with(info_snp, sqrt(2 * a1_af * (1 - a1_af)))
sd_ss <- with(info_snp, sqrt(2 * a1_sumdata_af * (1 - a1_sumdata_af)))
is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
df_beta <- info_snp[!is_bad, ]

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/GLGC/",trait), recursive=T)
setwd(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/GLGC/",trait))
tmp <- tempfile(tmpdir= "tmp-data")
for (chr in 1:22) {

  cat(chr, ".. ", sep = "")

  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

  corr_chr <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/LDpred2_format_LD/ref_geno/",race,"/LD_ref_chr", chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == 1) {
    ld <- Matrix::colSums(corr_chr^2)
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr_chr^2))
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = NCORES)))
ldsc_h2_est <- ldsc[["h2"]]
  h2_seq <- round(ldsc_h2_est * c(0.3, 0.7, 1, 1.4), 4)
  p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

dir.create("./ldpred2")
dir.create("./ldpred2/coef")
dir.create("./ldpred2/prs")
dir.create("./ldpred2/prs/tuning")
set.seed(1)  # to get the same result every time
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
save(beta_grid, df_beta, file=paste0("./ldpred2/coef/",race,"_ldpred2.RData"))

beta_grid[is.na(beta_grid)] = 0
beta_grid = data.frame(df_beta[,c('rsid','a1')], beta_grid)
colnames(beta_grid) = c('rsid','a1', paste0('e',1:nrow(params)))

for(chr in 1:22){
  print(chr)

  m <- df_beta$chr==chr
  prs.file <- beta_grid[m,,drop=F]
  write.table(prs.file, file = paste0("./ldpred2/coef/",race,"_ldpred2-chr",chr,".txt"),
              col.names = F,row.names = F,quote=F)

}



