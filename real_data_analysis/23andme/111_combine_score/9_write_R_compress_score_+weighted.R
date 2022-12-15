
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

########################################################
## 9_compress_score_+weighted.R

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/23andme/111_combine_score/9_compress_score_+weighted"))

a <- character()
b <- character()
n <- 0

#for (M in c(2,5)){
  for (targetethnic in c("EUR","AFR","AMR","EAS","SAS")){
    for (trait in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){

      n <- n+1
a <- paste0(a, n,"-",trait,"_",targetethnic,".Rout
")

b <- paste0(b,"trait='",trait,"' targetethnic='",targetethnic,"'
")
    }
  }
#}
writeLines(a,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/23andme/111_combine_score/9_compress_score_+weighted/ind_Rout.txt"))
writeLines(b,  paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/23andme/111_combine_score/9_compress_score_+weighted/ind_para.txt"))


a <- paste0("#!/usr/bin/env bash
#$ -N compress_score_from_single
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -m e
#$ -t 1-",n,"
#$ -M jzhan218@jhu.edu

readarray -t a < ind_Rout.txt
readarray -t b < ind_para.txt

module load conda_R/4.0

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../9_compress_score_+weighted.R ${a[$(($SGE_TASK_ID-1))]}
}
runr \"--args ${b[$(($SGE_TASK_ID-1))]}\"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/23andme/111_combine_score/9_compress_score_+weighted/ALL.sh"))


#qsub -l h=${skipnode} -hold_jid "cs23-*" ALL.sh
