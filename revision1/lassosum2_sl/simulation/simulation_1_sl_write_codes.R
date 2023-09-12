

setting=1

########################################################
## ../revision1/lassosum2_sl/simulation


dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/lassosum2_sl/simulation/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/lassosum2_sl/simulation/Setting_",setting))
#dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/"))


for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  #dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic))
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/lassosum2_sl/simulation/Setting_",setting,"/",ethnic))

  a <- character()
  b <- character()
  n <- 0

  for (rep in 1:10){
    #dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep))

  for (rho in 1:3){
    #for (size in 1:4){
    for (size in 1){
      n <- n+1
a <- paste0(a, n,"-rho_",rho,"_size_",size,"_rep_",rep,".Rout
")
b <- paste0(b, "ethnic='",ethnic,"' setting='",setting,"' rho='",rho,"' size='",size,"' rep='",rep,"'
")
}
  }
  }
  writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/lassosum2_sl/simulation/Setting_",setting,"/",ethnic,"/ind_Rout.txt"))
  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/lassosum2_sl/simulation/Setting_",setting,"/",ethnic,"/ind_para.txt"))

a <- paste0("#!/usr/bin/env bash
#$ -N four_tv_",ethnic,"_GA_",setting,"
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=1000G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../simulation_1_sl.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")
  writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/lassosum2_sl/simulation/Setting_",setting,"/",ethnic,"/simulation_1_sl.sh"))
}
