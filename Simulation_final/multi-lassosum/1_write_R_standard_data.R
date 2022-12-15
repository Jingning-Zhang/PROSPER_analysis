
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

########################################################

#ethnic=c("AFR","AMR","EAS","EUR","SAS")
ethnic2=paste0('c("',paste(ethnic,collapse = '","'),'")')
ethnic1=paste(ethnic,collapse = "_")

M <- length(ethnic); M

if(M==5){ memory <- 30 }else if(M==4){ memory <- 25 }else if(M==3){ memory <- 20 }else if(M==2){ memory <- 16 }

memory

#setting=2
#rep=1
#
#L=3
#Lc=3

########################################################
## 1_standard_data.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/Setting_",setting))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/1_standard_data"))

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data"))

a <- character()
b <- character()
n <- 0

for (rho in 1:3){
  for (size in 1:4){
    for (chr in 1:22){
         n <- n+1
a <- paste0(a, n,"-",M,"ethnicity-",ethnic1, "-rho",rho,"_size",size,"_rep_",rep,"_chr",chr,".Rout
")

b <- paste0(b,"rho='",rho,"' size='",size,"' chr='",chr,"' setting='",setting,"' rep='",rep,"' ethnic=",ethnic2,"
")
    }
  }
}
writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/1_standard_data/ind_Rout_rep_",rep,".txt"))
writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/1_standard_data/ind_para_rep_",rep,".txt"))


a <- paste0("#!/usr/bin/env bash
#$ -N five_mdata-",M,"ethnicity-",ethnic1,"-GA_",setting,"_rep_",rep,"
#$ -cwd
#$ -l mem_free=",memory+20,"G,h_vmem=",memory+20,"G,h_fsize=100G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout_rep_",rep,".txt
readarray -t b < ind_para_rep_",rep,".txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../../1_standard_data.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/1_standard_data/rep_",rep,".sh"))

