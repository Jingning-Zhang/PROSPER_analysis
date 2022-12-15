
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/1_LD_data/")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/1_LD_data/submitjob")


dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/LDdata/"))
dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/LDdata/",ethnic))

b <- paste0("#!/usr/bin/env bash
#$ -N multi-lassosum
#$ -cwd
#$ -m e
#$ -M jzhan218@jhu.edu
")

for (ethnic in c("AFR","AMR","EAS","EUR","SAS")){

  dir.create(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/LDdata/",ethnic))


a <- paste0("#!/usr/bin/env bash
#$ -N ld-",ethnic,"
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

Rscript ../1_compute_LD_by_chr.R \\
--packagedir /dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir \\
--refgeno /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/ref_chr \\
--snplist /dcs04/nilanjan/data/jzhang2/MEPRS/LD_simulation_GA/",ethnic,"/sumdata_snps_list.txt \\
--hg 19 \\
--chr $SGE_TASK_ID \\
--workdir /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/LDdata/",ethnic,"/


Rscript ../2_reformat_to_RData_by_chr.R \\
--packagedir /dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir \\
--hg 19 \\
--chr $SGE_TASK_ID \\
--workdir0 /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/LDdata/",ethnic,"/ \\
--workdir /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/LDdata/",ethnic,"/

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/1_LD_data/submitjob/",ethnic,".sh"))

  b <- paste0(b,"
qsub ",ethnic,".sh")

}

writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/multi-lassosum/1_LD_data/submitjob/ALL.sh"))


