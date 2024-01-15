
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/")

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge/0_compute_LD_UKBB//")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge/0_compute_LD_UKBB/1_standard_data/")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge/0_compute_LD_UKBB/1_standard_data/submitjob")

b <- paste0("#!/bin/bash
#SBATCH --job-name submit
#SBATCH --mem=1G
#SBATCH --exclude=compute-070,compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

")

for (ethnic in c("AFR","EAS","SAS")){

a <- paste0("#!/bin/bash
#SBATCH --job-name mdata-",ethnic,"
#SBATCH --mem=20G
#SBATCH --array=1-22
#SBATCH --exclude=compute-070,compute-089,compute-113,compute-116,compute-145,compute-147,compute-148
#SBATCH --mail-type=ALL --mail-user=jzhan218@jhu.edu

module load R

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge/0_compute_LD_UKBB/1_standard_data

mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block_UKBB/",ethnic,"/

Rscript ../1_compute_LD_by_chr.R \\
--packagedir /dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir \\
--refgeno /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/",ethnic,"/ref_chr \\
--hg 19 \\
--chr $SLURM_ARRAY_TASK_ID \\
--workdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block_UKBB/",ethnic,"/

Rscript ../2_reformat_to_RData_by_chr.R \\
--packagedir /dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir \\
--refgeno /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno_UKBB/",ethnic,"/ref_chr \\
--hg 19 \\
--chr $SLURM_ARRAY_TASK_ID \\
--workdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block_UKBB/",ethnic,"/

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge/0_compute_LD_UKBB/1_standard_data/submitjob/",ethnic,".sh"))

  b <- paste0(b,"
qsub -l h=${skipnode} ",ethnic,".sh")

}
writeLines(b, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/pacakge/0_compute_LD_UKBB/1_standard_data/submitjob/ALL.sh"))




