
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

########################################################
ethnic

ethnic2=paste0('c("',paste(ethnic,collapse = '","'),'")')
ethnic1=paste(ethnic,collapse = "_")

M <- length(ethnic)
M

########################################################
## 8_multi-lassosum_clean_score.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/8_multi-lassosum_clean_score_",para))

fixEUR='T'

#for (fixEUR in c("T","F")){

  if(fixEUR=="T"){
    path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/4_multi-lassosum_by_chr_",para,"_fixEUR")
  }else{
    path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para)
  }
  dir.create(paste0(path,"/Results/"))
  dir.create(paste0(path,"/Results/clean_results/"))


  if(fixEUR){ indfixEUR <- "_fixEUR" }else{ indfixEUR <- "" }

  a <- character()
  b <- character()
  n <- 0

  for (rho in 1:3){
    for (size in 1:4){
      n <- n+1
  a <- paste0(a, n,"-",M,"ethnicity-",ethnic1, "-rho",rho,"_size",size,"_rep_",rep,indfixEUR,".Rout
")

  b <- paste0(b,"rho='",rho,"' size='",size,"' setting='",setting,"' rep='",rep,"' ethnic=",ethnic2," L='",L,"' Lc='",Lc,"' para='",para,"' fixEUR='",fixEUR,"'
")
    }
  }

  writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/8_multi-lassosum_clean_score_",para,"/ind_Rout_rep_",rep,indfixEUR,".txt"))
  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/8_multi-lassosum_clean_score_",para,"/ind_para_rep_",rep,indfixEUR,".txt"))


  a <- paste0("#!/usr/bin/env bash
#$ -N seven_cs-",M,"ethnicity-",ethnic1,"-GA_",setting,"_rep_",rep,"_sub",indfixEUR,"
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout_rep_",rep,indfixEUR,".txt
readarray -t b <  ind_para_rep_",rep,indfixEUR,".txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../../8_multi-lassosum_clean_score.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")

  writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/8_multi-lassosum_clean_score_",para,"/rep_",rep,indfixEUR,".sh"))
#}

