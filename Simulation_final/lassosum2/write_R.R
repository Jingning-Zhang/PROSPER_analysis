
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

########################################################
## 1_prepare_LD.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/1_prepare_LD/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr"))


b <- paste0("#!/usr/bin/env bash
#$ -N one_LD
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){

a <- paste0("#!/usr/bin/env bash
#$ -N one_LD_",ethnic,"
#$ -cwd
#$ -pe local 8
#$ -l mem_free=20G,h_vmem=20G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../1_prepare_LD.R ",ethnic,".Rout
}
runr \"--args ethnic='",ethnic,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/1_prepare_LD/",ethnic,".sh"))

    b <- paste0(b,"
qsub ",ethnic,".sh")

}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/1_prepare_LD/ALL.sh"))




########################################################
## 1_prepare_LD_add_together.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/1_prepare_LD_add_together/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr"))


b <- paste0("#!/usr/bin/env bash
#$ -N one_LD
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){

a <- paste0("#!/usr/bin/env bash
#$ -N one_LD_",ethnic,"
#$ -cwd
#$ -pe local 6
#$ -l mem_free=20G,h_vmem=20G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../1_prepare_LD_add_together.R ",ethnic,".Rout
}
runr \"--args ethnic='",ethnic,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/1_prepare_LD_add_together/",ethnic,".sh"))

    b <- paste0(b,"
qsub ",ethnic,".sh")

}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/1_prepare_LD_add_together/ALL.sh"))




########################################################
## 2_prepare_sumdata.R

setting=5

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/2_prepare_sumdata"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr"))

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",ethnic))
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/2_prepare_sumdata/",ethnic))

  a <- character()
  b <- character()
  n <- 0

  for (rep in 1:10){
  for (rho in 1:3){
    for (size in 1:4){
      n <- n+1
a <- paste0(a, n,"-rho_",rho,"_size_",size,"_rep_",rep,".Rout
")
b <- paste0(b, "ethnic='",ethnic,"' setting='",setting,"' rho='",rho,"' size='",size,"' rep='",rep,"'
")
    }
  }
  }
  writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/2_prepare_sumdata/",ethnic,"/ind_Rout.txt"))
  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/2_prepare_sumdata/",ethnic,"/ind_para.txt"))

a <- paste0("#!/usr/bin/env bash
#$ -N two_sumdata_",ethnic,"_GA_",setting,"
#$ -cwd
#$ -pe local 2
#$ -l mem_free=10G,h_vmem=10G,h_fsize=1000G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../../2_prepare_sumdata.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")
  writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/2_prepare_sumdata/",ethnic,"/2_prepare_sumdata.sh"))
}

########################################################
## 3_run_lassosum2.R

setting=5

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/3_run_lassosum2"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/"))

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic))

  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/3_run_lassosum2/",ethnic))

  a <- character()
  b <- character()
  n <- 0

  for (rep in 1:10){
      dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep))

  for (rho in 1:3){
    for (size in 1:4){
      n <- n+1
a <- paste0(a, n,"-rho_",rho,"_size_",size,"_rep_",rep,".Rout
")
b <- paste0(b, "ethnic='",ethnic,"' setting='",setting,"' rho='",rho,"' size='",size,"' rep='",rep,"'
")
}
  }
  }
  writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/3_run_lassosum2/",ethnic,"/ind_Rout.txt"))
  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/3_run_lassosum2/",ethnic,"/ind_para.txt"))

a <- paste0("#!/usr/bin/env bash
#$ -N three_run_",ethnic,"_GA_",setting,"
#$ -cwd
#$ -pe local 6
#$ -l mem_free=15G,h_vmem=15G,h_fsize=1000G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../../3_run_lassosum2.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")
  writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/3_run_lassosum2/",ethnic,"/3_run_lassosum2.sh"))
}



########################################################
## 4_tuning_and_validation.R


dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/"))


for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic))
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation/",ethnic))

  a <- character()
  b <- character()
  n <- 0

  for (rep in 1:10){
    dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep))

  for (rho in 1:3){
    for (size in 1:4){
      n <- n+1
a <- paste0(a, n,"-rho_",rho,"_size_",size,"_rep_",rep,".Rout
")
b <- paste0(b, "ethnic='",ethnic,"' setting='",setting,"' rho='",rho,"' size='",size,"' rep='",rep,"'
")
}
  }
  }
  writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation/",ethnic,"/ind_Rout.txt"))
  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation/",ethnic,"/ind_para.txt"))

a <- paste0("#!/usr/bin/env bash
#$ -N four_tv_",ethnic,"_GA_",setting,"
#$ -cwd
#$ -pe local 12
#$ -l mem_free=2G,h_vmem=2G,h_fsize=1000G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../../4_tuning_and_validation.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")
  writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation/",ethnic,"/4_tuning_and_validation.sh"))
}


########################################################
## 4_tuning_and_validation_EUR.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation_EUR"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/"))


for (ethnic in c("AFR","AMR","EAS","SAS")){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic))

  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation_EUR/",ethnic))

  a <- character()
  b <- character()
  n <- 0

  for (rep in 1:10){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep))

  for (rho in 1:3){
      n <- n+1
a <- paste0(a, n,"-rho_",rho,"_rep_",rep,".Rout
")
b <- paste0(b, "ethnic='",ethnic,"' setting='",setting,"' rho='",rho,"' rep='",rep,"'
")
  }
  }
  writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation_EUR/",ethnic,"/ind_Rout.txt"))
  writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation_EUR/",ethnic,"/ind_para.txt"))

a <- paste0("#!/usr/bin/env bash
#$ -N four_tvEUR_",ethnic,"_GA_",setting,"
#$ -cwd
#$ -pe local 12
#$ -l mem_free=2G,h_vmem=2G,h_fsize=1000G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../../4_tuning_and_validation_EUR.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")
  writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/4_tuning_and_validation_EUR/",ethnic,"/4_tuning_and_validation_EUR.sh"))
}


########################################################
## 5_weighted_lassosum2.R


dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/5_weighted_lassosum2"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/"))

a <- character()
b <- character()
n <- 0

for (ethnic in c("AFR","AMR","EAS","SAS")){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic))

  for(rep in 1:10){

  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/Setting_",setting,"_lassosum2/",ethnic,"/rep_",rep))
  n <- n+1
a <- paste0(a, n,"-",ethnic,"_rep_",rep,".Rout
")
b <- paste0(b, "ethnic='",ethnic,"' setting='",setting,"' rep='",rep,"'
")

}
}
writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/5_weighted_lassosum2/ind_Rout.txt"))
writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/5_weighted_lassosum2/ind_para.txt"))

a <- paste0("#!/usr/bin/env bash
#$ -N five_weightedlassosum2_GA_",setting,"
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=1000G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../5_weighted_lassosum2.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/lassosum2/Setting_",setting,"/5_weighted_lassosum2/5_weighted_lassosum2.sh"))
