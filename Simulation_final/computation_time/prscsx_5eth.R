
args = commandArgs(trailingOnly = T)
repl = as.numeric(args[[1]])
print(repl)


start_time <- Sys.time()

library(bigreadr)
library(stringr)
library(dplyr)


path_to_sum = "/dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prscsx/sumdata/"
out.dir = paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prscsx/results/",repl)

dir.create(out.dir)

phi = c(1,1E-02,1E-04,1E-06)

set.seed(repl)

for(k in 1:length(phi)){

  print(k)


  system(paste0("export MKL_NUM_THREADS=1; export NUMEXPR_NUM_THREADS=1; export OMP_NUM_THREADS=1;
                python /dcs04/nilanjan/data/jzhang2/TOOLS/prscsx/PRScsx/PRScsx.py",
              " --ref_dir=/dcs04/nilanjan/data/jzhang2/TOOLS/prscsx/LDref/",
              " --bim_prefix=",path_to_sum,"alleth",
              " --sst_file=",path_to_sum,"EUR.txt,",path_to_sum,"AFR.txt,",path_to_sum,"AMR.txt,",path_to_sum,"EAS.txt,",path_to_sum,"SAS.txt",
              " --n_gwas=100000,15000,15000,15000,15000",
              " --pop=EUR,AFR,AMR,EAS,SAS",
              " --chrom=22",
              " --phi=",phi[k],
              " --out_dir=",out.dir,
              " --out_name=five"))

}

####################################
## clean score by tuning parameter settings

ethnic <- c("EUR","AFR","AMR","EAS","SAS")
for ( j in 1:length(ethnic)){
  eth <- ethnic[j]
  tmp <- dir(out.dir)
  tmp <- tmp[str_detect(tmp, paste0("five_",eth)) & (str_detect(tmp,"chr22.txt"))]
  for(k in 1:length(tmp)){
    dat <- fread2(paste0(out.dir, "/",tmp[k]))[,c(2,4,5,6)]
    if(k==1) { score <- dat } else score <- full_join(score, dat, by=c("V2","V4","V5"))
  }
  fwrite2(score, paste0(out.dir,"/five_",eth,"_alltuning.txt"), col.names = F, sep="\t", nThread=1)
  n_col <- ncol(score)

  system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 ",
                "--score-col-nums 4-",n_col," --threads 1 ",
                "--score ",out.dir,"/five_",eth,"_alltuning.txt cols=+scoresums,-scoreavgs ",
                "--bfile /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno ",
                "--out ",out.dir,"/five_",eth,"_alltuning"))



## match up phenotype

  fam <- read.table(paste("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno.fam",sep=''),as.is=T)
  pheno <- read.table("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/pheno.fam", as.is=T)
  m <- match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
  m.keep <- !is.na(m)
  fam <- fam[m.keep,,drop=F]
  pheno <- pheno[m[m.keep],,drop=F]
  SCORE <- fread2(paste0(out.dir,"/five_",eth,"_alltuning.sscore"))
  m <- match( paste(fam[,1],fam[,2]) , paste(SCORE[,1],SCORE[,2]) )
  m.keep <- !is.na(m)
  fam <- fam[m.keep,]
  pheno <- pheno[m.keep,]
  SCORE <- SCORE[m[m.keep],]
  SCORE_id <- SCORE[,1:4]
  SCORE <- SCORE[,-1:-4]
  colnames(SCORE) <- paste0("score",1:ncol(SCORE))
  R2 <- numeric(length = ncol(SCORE))
  for (i in 1:ncol(SCORE)){
    fit <- lm( pheno[,3] ~ SCORE[,i] )
    R2[i] <- summary(fit)$r.square
  }
  best <- which.max(R2)
  best_score_tmp <- score[,c(1,2,3,best+3)]
  colnames(best_score_tmp)[4] <- "weight"
  best_score_tmp <- best_score_tmp[best_score_tmp$weight!=0,]
  fwrite2(best_score_tmp, paste0(out.dir,"/five_",eth,"_best_score.txt"), col.names = T, sep="\t", nThread=1)

  res_tmp <- cbind(SCORE_id,SCORE[,best])
  colnames(res_tmp)[5] <- eth

  if(j==1){
    res <- res_tmp[,c(1,2,5)]
    best_score <- best_score_tmp
  }else{
    res <- inner_join(res, res_tmp[,c(1,2,5)])
    best_score <- full_join(best_score, best_score_tmp, by=c("V2","V4","V5"))
  }
}
colnames(best_score) <- c("rsid","a1","a0",ethnic)
best_score[is.na(best_score)] <- 0


####################################
## weighted by ancestry

m <- match( paste(fam[,1],fam[,2]) , paste(res[,1],res[,2]) )
m.keep <- !is.na(m)
fam <- fam[m.keep,]
pheno <- pheno[m.keep,]
res <- res[m[m.keep],]
SCORE_id <- res[,1:2]
SCORE <- res[,-1:-2]
fit <- lm( pheno ~. , data=cbind(pheno=pheno[,3],SCORE ))

tmp <- rbind(rep(0,ncol(SCORE)), diag(ncol(SCORE)))
colnames(tmp) <- colnames(data.frame(SCORE))
tmp <- predict(fit, as.data.frame(tmp))
coef <- (tmp-tmp[1])[-1]

####################################
## summarize final scores

ensemble_score <- data.frame(best_score[,1:3], weight=as.matrix(best_score[,-(1:3)]) %*% matrix(coef, ncol=1))
ensemble_score <- ensemble_score[ensemble_score$weight!=0,]

fwrite2(ensemble_score, paste0(out.dir,"/five_final_score.txt"), col.names = T, sep="\t", nThread=1)

end_time <- Sys.time()

runtime <- as.numeric(end_time-start_time, units="secs")

saveRDS(runtime, file= paste0(out.dir,"/five_runtime.rds"))


