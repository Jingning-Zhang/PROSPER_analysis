
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }


ethnic2=paste0("c('",paste(ethnic,collapse = "','"),"')")
ethnic1=paste(ethnic,collapse = "_")

M <- length(ethnic)

if(M==5){
  memory <- 20
}else if(M==4){
  memory <- 16
}else if(M==3){
  memory <- 15
}else if(M==2){
  memory <- 10
}

########################################################


dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum/",ethnic1))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum/",ethnic1,"/1_standard_data"))

dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/multi-lassosum/",ethnic1))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/multi-lassosum/",ethnic1,"/1_standard_data"))

########################################################
## 1_standard_data

b <- paste0("#!/usr/bin/env bash
#$ -N multi-lassosum",ethnic1,"
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

")

#for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
for (trait in c("HDL")){
  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/multi-lassosum/",ethnic1,"/1_standard_data/",trait))

a <- paste0("#!/usr/bin/env bash
#$ -N mdata-",M,"ethnicity-",ethnic1,"_",trait,"
#$ -cwd
#$ -l mem_free=",memory,"G,h_vmem=",memory,"G,h_fsize=100G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

module load R/3.6.1

runr(){
    R CMD BATCH --no-save --no-restore \"$1\"  ../../1_standard_data.R ",ethnic1,"_",trait,"_chr$SGE_TASK_ID.Rout
}
runr \"--args ethnic=",ethnic2," trait='",trait,"' chr='$SGE_TASK_ID' \"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum/",ethnic1,"/1_standard_data/",trait,".sh"))

  b <- paste0(b,"
qsub -l h=${skipnode} -hold_jid 'mdata-*_",trait,"' ",trait,".sh")

}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum/",ethnic1,"/1_standard_data/ALL.sh"))

