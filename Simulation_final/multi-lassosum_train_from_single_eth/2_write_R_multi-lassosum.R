
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

d <- floor(L*Lc/6000)+1
if(M==5){ memory <- 20*d }else if(M==4){ memory <- 17*d }else if(M==3){ memory <- 13*d }else if(M==2){ memory <- 10*d }

########################################################
## 2_MultiEthnic_runPRS_by_chr.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/2_multi-lassosum_by_chr_",para))

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para,"/tuning"))


a <- character()
b <- character()
n <- 0

for (rho in 1:3){
  for (size in 1:4){
    for (chr in 1:22){
         n <- n+1
a <- paste0(a, n,"-",M,"ethnicity-",ethnic1, "-rho",rho,"_size",size,"_rep_",rep,"_chr",chr,".Rout
")

b <- paste0(b,"rho='",rho,"' size='",size,"' chr='",chr,"' setting='",setting,"' rep='",rep,"' ethnic=",ethnic2," L='",L,"' Lc='",Lc,"' para='",para,"'
")
    }
  }
}
writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/ind_Rout_rep_",rep,".txt"))
writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/ind_para_rep_",rep,".txt"))


a <- paste0("#!/usr/bin/env bash
#$ -N six_meprs-",M,"ethnicity-",ethnic1,"-GA_",setting,"_rep_",rep,"_sub
#$ -cwd
#$ -l mem_free=",memory,"G,h_vmem=",memory,"G,h_fsize=100G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout_rep_",rep,".txt
readarray -t b < ind_para_rep_",rep,".txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../../2_multi-lassosum.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/2_multi-lassosum_by_chr_",para,"/rep_",rep,".sh"))


