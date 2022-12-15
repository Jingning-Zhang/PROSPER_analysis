
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum/")

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum/0_compute_LD//")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum/0_compute_LD/1_standard_data/")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum/0_compute_LD/1_standard_data/submitjob")

b <- paste0("#!/usr/bin/env bash
#$ -N multi-lassosum
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-048|compute-053|compute-054|compute-066|compute-068|compute-070|compute-076|compute-113)'

")

for (ethnic in c("AFR","AMR","EAS","EUR","SAS")){
     dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum/",ethnic))

for (trait in c("HDL","LDL","logTG","nonHDL","TC")){

  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum/",ethnic,"/",trait))

a <- paste0("#!/usr/bin/env bash
#$ -N mdata-",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum/0_compute_LD/1_standard_data

Rscript ../1_compute_LD_by_chr.R \\
--packagedir /dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir \\
--refgeno /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/",ethnic,"/",trait,"/ref_chr \\
--hg 19 \\
--chr $SGE_TASK_ID \\
--workdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum/",ethnic,"/",trait,"


Rscript ../2_reformat_to_RData_by_chr.R \\
--packagedir /dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir \\
--refgeno /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/",ethnic,"/",trait,"/ref_chr \\
--hg 19 \\
--chr $SGE_TASK_ID \\
--workdir /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/multi-lassosum/",ethnic,"/",trait,"

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum/0_compute_LD/1_standard_data/submitjob/",ethnic,"_",trait,".sh"))

  b <- paste0(b,"
qsub -l h=${skipnode} ",ethnic,"_",trait,".sh")

}
}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/multi-lassosum/0_compute_LD/1_standard_data/submitjob/ALL.sh"))



#for (ethnic in c("AFR","AMR","EAS","EUR","SAS")){
#
#for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
#  a <- paste0("mv /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/final_results/",ethnic,"_",trait,"_validation_R2_best.rds /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/lassosum2/Results/final_results/",ethnic,"_",trait,"_best.rds")
#}}
