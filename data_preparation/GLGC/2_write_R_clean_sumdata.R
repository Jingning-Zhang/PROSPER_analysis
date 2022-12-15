
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

memory <- 20

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/data_preparation/GLGC/2_clean_sumdata"))

a <- character()
b <- character()
n <- 0
for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
  for (race in c("AFR","AMR","EAS","EUR","SAS")){
         n <- n+1
a <- paste0(a, n,"-",trait,"-",race,".Rout
")

b <- paste0(b,"trait='",trait,"' race='",race,"'
")

  }
}
writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/data_preparation/GLGC/2_clean_sumdata/ind_Rout.txt"))
writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/data_preparation/GLGC/2_clean_sumdata/ind_para.txt"))


a <- paste0("#!/usr/bin/env bash
#$ -N clean_sumdata_glgc
#$ -cwd
#$ -l mem_free=",memory,"G,h_vmem=",memory,"G,h_fsize=100G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../2_clean_sumdata.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/data_preparation/GLGC/2_clean_sumdata/submit.sh"))


