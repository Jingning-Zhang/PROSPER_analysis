
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

########################################################
## 1_prepare_data_lassosum2.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/1_prepare_data_lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/"))



b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){


a <- paste0("#!/usr/bin/env bash
#$ -N data_",ethnic,"_",trait,"
#$ -cwd
#$ -pe local 8
#$ -l mem_free=20G,h_vmem=20G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../1_prepare_data_lassosum2.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/1_prepare_data_lassosum2/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/1_prepare_data_lassosum2/ALL.sh"))




########################################################
## 2_run_lassosum2.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/2_run_lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/"))


b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){


a <- paste0("#!/usr/bin/env bash
#$ -N run_",ethnic,"_",trait,"
#$ -cwd
#$ -pe local 8
#$ -l mem_free=20G,h_vmem=20G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../2_run_lassosum2.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/2_run_lassosum2/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/2_run_lassosum2/ALL.sh"))


########################################################
## 3_lassosum2_clean_score.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/3_lassosum2_clean_score/"))
dir.create(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/"))


b <- paste0("#!/usr/bin/env bash
#$ -N lassosum2
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){


a <- paste0("#!/usr/bin/env bash
#$ -N run_",ethnic,"_",trait,"
#$ -cwd
#$ -pe local 2
#$ -l mem_free=15G,h_vmem=15G,h_fsize=1000G
#$ -m e
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../3_lassosum2_clean_score.R ",ethnic,"_",trait,".Rout
}
runr \"--args ethnic='",ethnic,"' trait='",trait,"'\"

")
writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/3_lassosum2_clean_score/",ethnic,"_",trait,".sh"))

    b <- paste0(b,"
qsub ",ethnic,"_",trait,".sh")

  }
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/lassosum2/3_lassosum2_clean_score/ALL.sh"))

