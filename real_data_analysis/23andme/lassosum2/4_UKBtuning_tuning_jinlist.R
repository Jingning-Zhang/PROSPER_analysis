
#library(readr)
#p <- read_csv("/dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/EUR/pheno.txt")
#set.seed(1)
#a <- sort(sample(1:nrow(p),floor(nrow(p)/2)))
#b <- (1:nrow(p))[-a]
#
#write.table(data.frame(p$eid,p$eid),
#            quote=F, row.names=F, col.names=F,
#            "/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning_jinlist/allsample.eid")
#
#write.table(data.frame(p$eid[a],p$eid[a]),
#            quote=F, row.names=F, col.names=F,
#            "/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning_jinlist/tuning.eid")
#
#write.table(data.frame(p$eid[b],p$eid[b]),
#            quote=F, row.names=F, col.names=F,
#            "/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/UKBtuning_jinlist/validation.eid")

############################################################
############################################################
## lassosum2

library(readr)

system(paste0("awk '{print $1}' /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum2.score",
              " > /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum2.snp"))

tmp <- paste0("#!/usr/bin/env bash
#$ -N lassosum2_EUR
#$ -cwd
#$ -pe local 8
#$ -l mem_free=1G,h_vmem=1G,h_fsize=1000G
#$ -t 1-22
#$ -m e
#$ -M jzhan218@jhu.edu

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 8 --rm-dup exclude-all",
            " --bfile /dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/EUR/chr${SGE_TASK_ID}",
            " --extract /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum2.snp",
            " --score /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum2.score",
            " cols=+scoresums,-scoreavgs --score-col-nums 3-202",
            " --out /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score_lassosum2/EUR_height_chr${SGE_TASK_ID}"
)

writeLines(tmp,"/dcs04/nilanjan/data/23andme/Analysis/JZ/codesR/UKB_EUR_height_lassosum2.sh")


### lassosum
#
#system(paste0("awk '{print $1}' /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum.score",
#              " > /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum.snp"))
#
#tmp <- paste0("#!/usr/bin/env bash
##$ -N lassosum
##$ -cwd
##$ -pe local 15
##$ -l mem_free=1G,h_vmem=1G,h_fsize=1000G
##$ -t 1-22
##$ -m e
##$ -M jzhan218@jhu.edu
#
#/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 15 --rm-dup exclude-all",
#            " --bfile /dcs04/nilanjan/data/jjin/UKB/ldpred_validation_mega/EUR/chr${SGE_TASK_ID}",
#            " --extract /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum.snp",
#            " --score /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score/EUR_height_lassosum.score",
#            " cols=+scoresums,-scoreavgs --score-col-nums 3-202",
#            " --out /dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score_lassosum/EUR_height_chr${SGE_TASK_ID}"
#)
#
#writeLines(tmp,"/dcs04/nilanjan/data/23andme/Analysis/JZ/codesR/UKB_EUR_height_lassosum.sh")

############################################################
############################################################

library(readr)
for (chr in 1:22){
  print(chr)
  score <- read_tsv(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score_lassosum2/EUR_height_chr",chr,".sscore"), col_types = cols())
  if(chr==1){
    ID <- score[,1:4]
    SCORE <- score[, -1:-4]
  }else{
    SCORE <- SCORE + score[, -1:-4]
  }
}
SCORE <- cbind(ID, SCORE)
write_tsv(SCORE, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score_lassosum2/EUR_height_allchr.sscore"))


#for (chr in 1:22){
#  print(chr)
#  score <- read_tsv(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score_lassosum/EUR_height_chr",chr,".sscore"), col_types = cols())
#  if(chr==1){
#    ID <- score[,1:4]
#    SCORE <- score[, -1:-4]
#  }else{
#    SCORE <- SCORE + score[, -1:-4]
#  }
#}
#SCORE <- cbind(ID, SCORE)
#write_tsv(SCORE, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/UKBtuning_jinlist/score_lassosum/EUR_height_allchr.sscore"))





