
rm(list=ls())
setwd('~/prs/')
library(data.table)
library(Rcpp)
library(readr)
library(MASS) # for mvrnorm
library(reshape) # for melt
library(parallel)
library(RcppTN)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(genio) # for read_plink
library(dplyr)
library(stringr)
library(gdata)
library(R.utils) # for gzip
library(pROC)
library(bigsnpr)
library(bigparallelr)

races = c('EUR','AFR','AMR')
traits = c('height','bmi')

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = races[as.numeric(temp1)]
trait = traits[as.numeric(temp2)]

ldr = 3/1000
ncores = 1

temdir = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait)
if (!dir.exists(temdir)){dir.create(temdir)}
temdir = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait,'/',race)
if (!dir.exists(temdir)){dir.create(temdir)}

tedir = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait,'/',race,'/intermediate/')
if (!dir.exists(tedir)){dir.create(tedir)}

temdir = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait,'/',race,'/ldpred2/')
if (!dir.exists(temdir)){dir.create(temdir)}

NCORES <-  nb_cores()


  temfile = paste0(temdir,"ldpred2-chr",chr,".txt")
  outfile = temfile
  sumraw = bigreadr::fread2(paste0("/dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/",race,"/",trait,".txt"))
  sumraw = sumraw[,c('rsID','CHR','POS_b38','N','BETA','SE','P','A1','A2')]
  colnames(sumraw) = c('SNP_ID','CHR','POS','N','BETA','SE','PVAL','REF','ALT') # REF: allele corresponding to BETA

  valdat = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr,'.bim'))
  sumraw = sumraw[sumraw$SNP_ID %in% valdat[,2],]

  #### read in reference data
  refdir = paste0('/dcs04/nilanjan/data/jjin/prs/realdata/glgc/refgeno/',race)
  #if (!dir.exists(refdir)){dir.create(refdir)}
  #temfile = paste0(refdir,'/chr',chr,'.bk')
  #system(paste0('rm -rf ',temfile))
  #temfile = paste0(refdir,'/chr',chr,'.rds')
  #system(paste0('rm -rf ',temfile))
  #snp_readBed(paste0(refdir,'/chr',chr,'.bed'))
  obj.bigSNP <- snp_attach(paste0(refdir,'/chr',chr,'.rds'))
  map <- obj.bigSNP$map[-c(3)]
  names(map) <- c("chr", "rsid", "pos", "a0", "a1") # a1 - alt # c("chr", "pos", "a0", "a1")

  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait,'/',race,'/intermediate/'), ncores = ncores)

  sumstats = sumraw[sumraw$CHR == chr,c('CHR', 'SNP_ID', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'N')]
  set.seed(2020)
  names(sumstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")
  sumstats = sumstats[sumstats$rsid %in% map$rsid,]
  info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F) # important: for real data, strand_flip = T
  rownames(info_snp) = info_snp$rsid
  ind.chr <- which(info_snp$chr == chr)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

  corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = 3,
                   infos.pos = POS2[ind.chr2], size =  ldr)
  corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))

  # Automatic model
  ldsc <- snp_ldsc2(corr0, df_beta)
  h2_est <- ldsc[["h2"]]
  print(paste0('Complete data preparation'))

  h2_seq <- c(0.3, 0.7, 1, 1.4)
  p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  params <- expand.grid(p = p_seq, h2 = signif(abs(h2_est) * h2_seq, 3), sparse = c(FALSE))

  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, burn_in = 50, num_iter = 200, ncores = NCORES)
  rownames(beta_grid) = info_snp$rsid
  beta_grid = cbind(beta_grid, info_snp[,c('a0','a1','rsid')])
  colnames(beta_grid) = c(paste0('e',1:nrow(params)), 'a0','a1','rsid')

  beta_grid[is.na(beta_grid)] = 0
  beta_grid = as.data.frame(beta_grid)
  write_delim(beta_grid,file = outfile, delim='\t')
  rm(corr0, corr)
  print(paste0('Complete'))



  tedir = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait,'/',race,'/ldpred2/effect/')
  if (!dir.exists(tedir)) dir.create(tedir)
  tedir = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait,'/',race,'/ldpred2/score/')
  if (!dir.exists(tedir)) dir.create(tedir)
  outdir = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait,'/',race,'/ldpred2/')

  outfile = paste0(outdir, 'score/ldpred2-chr', chr,'.sscore')

  # -------- PRS:
  prs.file = data.frame(SNP = beta_grid$rsid, ALT = beta_grid$a0, BETA = beta_grid[,1:(ncol(beta_grid)-3)])
  tem = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr,'.bim'))
  dupid = tem[duplicated(tem[,2]),2]
  prs.file = prs.file[!(prs.file$SNP %in% dupid),]
  write.table(prs.file,file = paste0(outdir,'effect/ldpred2-chr',chr,'.txt'),
              col.names = T,row.names = F,quote=F)
  prscode = paste(paste0('/dcl01/chatterj/data/jin/software/plink2'),
                  paste0('--score ',  outdir,'effect/ldpred2-chr',chr,'.txt'),
                  'cols=+scoresums,-scoreavgs',
                  paste0('--score-col-nums 3-',ncol(prs.file)),
                  paste0(' --bfile /dcs04/nilanjan/data/jjin/UKBval/',race,'/genotype/chr',chr),
                  " --threads 1",
                  paste0(' --out ', outdir, 'score/ldpred2-chr', chr))

  system(prscode)

  print(paste0('Completed ',trait,race,chr))
