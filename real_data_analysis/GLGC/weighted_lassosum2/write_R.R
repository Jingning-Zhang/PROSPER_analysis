
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

########################################################
## 2_test.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/weighted_lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/weighted_lassosum2/2_test/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/weighted_lassosum2/"))


b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){

    dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/weighted_lassosum2/test_all_prs/"))
    dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/weighted_lassosum2/final_results/"))

a <- paste0("#!/usr/bin/env bash
#$ -N weighted_lassosum2_",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../2_test.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/weighted_lassosum2/2_test/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub -l h=${skipnode} ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/weighted_lassosum2/2_test/ALL.sh"))

## -hold_jid clean_score