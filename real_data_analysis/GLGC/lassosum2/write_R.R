
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

########################################################
## 1_prepare_data_lassosum2.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/1_prepare_data_lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/"))



b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){


a <- paste0("#!/usr/bin/env bash
#$ -N data_",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=40G,h_vmem=40G,h_fsize=1000G
#$ -pe local 2
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../1_prepare_data_lassosum2.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/1_prepare_data_lassosum2/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub -l h=${skipnode} -hold_jid ",ethnic,"_",trait," ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/1_prepare_data_lassosum2/ALL.sh"))



########################################################
## 2_run_lassosum2.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/2_run_lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/"))


b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){


a <- paste0("#!/usr/bin/env bash
#$ -N run_",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=70G,h_vmem=70G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../2_run_lassosum2.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/2_run_lassosum2/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub -l h=${skipnode} -hold_jid data_",ethnic,"_",trait," ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/2_run_lassosum2/ALL.sh"))



########################################################
## 3_lassosum2_clean_score.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/3_lassosum2_clean_score/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/"))


b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){


a <- paste0("#!/usr/bin/env bash
#$ -N clean_",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../3_lassosum2_clean_score.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/3_lassosum2_clean_score/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub -l h=${skipnode} -hold_jid run_",ethnic,"_",trait," ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/3_lassosum2_clean_score/ALL.sh"))



########################################################
## 4_ukbb_tuning+validation.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/4_ukbb_tuning+validation/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/"))


b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){

a <- paste0("#!/usr/bin/env bash
#$ -N tv_",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../4_ukbb_tuning+validation.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/4_ukbb_tuning+validation/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub -l h=${skipnode} -hold_jid clean_",ethnic,"_",trait," ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/4_ukbb_tuning+validation/ALL.sh"))



########################################################
## 5_ukbb_tuning+validation_best_EUR.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/5_ukbb_tuning+validation_best_EUR/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/"))


b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){


a <- paste0("#!/usr/bin/env bash
#$ -N bestEUR_",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../5_ukbb_tuning+validation_best_EUR.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/5_ukbb_tuning+validation_best_EUR/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub -l h=${skipnode} -hold_jid tv_",ethnic,"_",trait," ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/5_ukbb_tuning+validation_best_EUR/ALL.sh"))








########################################################
## 1_prepare_data_lassosum2.R check err

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/1_prepare_data_lassosum2_check_err/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/"))


b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){


a <- paste0("#!/usr/bin/env bash
#$ -N data_",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=70G,h_vmem=70G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../1_prepare_data_lassosum2_check_err.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/1_prepare_data_lassosum2_check_err/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub -l h=${skipnode} -hold_jid ",ethnic,"_",trait," ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/1_prepare_data_lassosum2_check_err/ALL.sh"))


