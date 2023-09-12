
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }
L <- as.integer(L); Lc <- as.integer(Lc);

########################################################

ethnic2=paste0('c("',paste(ethnic,collapse = '","'),'")')
ethnic1=paste(ethnic,collapse = "_")

M <- length(ethnic)
M

#d <- floor(L^M*Lc/6000)+1
#if(M==5){ memory <- 50*d }else if(M==4){ memory <- 45*d }else if(M==3){ memory <- 40*d }else if(M==2){ memory <- 25*d }

memory <- 20

########################################################
## 2_MultiEthnic_runPRS_by_chr.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/",ethnic1))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para))

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/",ethnic1))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/tuning"))


a <- character()
b <- character()
n <- 0

#for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
for (trait in c("HDL")){

    for (chr in 1:22){
         n <- n+1
a <- paste0(a, n,"-",M,"ethnicity-",ethnic1, "-",trait,"_chr",chr,".Rout
")

b <- paste0(b,"trait='",trait,"' chr='",chr,"' ethnic=",ethnic2," L='",L,"' Lc='",Lc,"' para='",para,"'
")

  }
}
writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/ind_Rout.txt"))
writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/ind_para.txt"))


a <- paste0("#!/usr/bin/env bash
#$ -N meprs-",M,"ethnicity-",ethnic1,"
#$ -cwd
#$ -l mem_free=",memory,"G,h_vmem=",memory,"G,h_fsize=100G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../2_multi-lassosum.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/ALL.sh"))

