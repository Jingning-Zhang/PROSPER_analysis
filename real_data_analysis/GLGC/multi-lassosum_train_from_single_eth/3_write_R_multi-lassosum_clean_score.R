
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
#if(M==5){ memory <- 40*d }else if(M==4){ memory <- 37*d }else if(M==3){ memory <- 25*d }else if(M==2){ memory <- 30*d }

memory <- 20

########################################################
## 3_multi-lassosum_clean_score.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/3_multi-lassosum_by_chr_",para,"_clean_score"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/Results"))

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/prelim_results_table/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/prelim_results_table/",M,"ethnicity"))


a <- character()
b <- character()
n <- 0

for (trait in c("HDL","LDL","logTG","nonHDL","TC")){

  for (targetethnic in ethnic){
    #if(targetethnic == "EUR"){ next }
    dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/prelim_results_table/",M,"ethnicity/",targetethnic))
    dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/prelim_results_table/",M,"ethnicity/",targetethnic,"/",trait))

      n <- n+1
a <- paste0(a, n,"-",M,"ethnicity-",ethnic1, "-",trait,"_",targetethnic,".Rout
")

b <- paste0(b,"trait='",trait,"' ethnic=",ethnic2," L='",L,"' Lc='",Lc,"' para='",para,"' targetethnic='",targetethnic,"'
")

  }
}
writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/3_multi-lassosum_by_chr_",para,"_clean_score/ind_Rout.txt"))
writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/3_multi-lassosum_by_chr_",para,"_clean_score/ind_para.txt"))


a <- paste0("#!/usr/bin/env bash
#$ -N cs-",M,"ethnicity-",ethnic1,"
#$ -cwd
#$ -l mem_free=",memory,"G,h_vmem=",memory,"G,h_fsize=100G
#$ -pe local 1
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../3_multi-lassosum_clean_score.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/3_multi-lassosum_by_chr_",para,"_clean_score/ALL.sh"))
paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum_train_from_single_eth/",ethnic1,"/3_multi-lassosum_by_chr_",para,"_clean_score/ALL.sh")




