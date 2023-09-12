
################################################################
rm(list=ls())
args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }


a <- read.table("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/unrelated_whites.id",header=T)
b <- read.table("/dcs04/nilanjan/data/jzhang2/UKBB/ancestry_prediction/final/tuning+validation/EUR.all_id",header=T)

a <- a[!(a$FID %in% b$FID),]
set.seed(1)
b <- a[sample(1:nrow(a), 10000),]
readr::write_tsv(b, "/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/UKBB_EUR.id")


"mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/SAS/
cp -r /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/SAS/HDL /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/SAS/
"

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum/0_compute_LD/0_EUR_UKBB")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/EUR")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/EUR/HDL")

a <- paste0("#!/usr/bin/env bash
#$ -N ukbb_ld_eur
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=100G
#$ -m e
#$ -t 1-22
#$ -M jzhan218@jhu.edu

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 \\
--pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr$SGE_TASK_ID \\
--extract /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/snps/EUR_HDL_rsid.txt \\
--keep /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/UKBB_EUR.id \\
--make-bed \\
--rm-dup exclude-all \\
--threads 1 \\
--out /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/EUR/HDL/ref_chr$SGE_TASK_ID

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum/0_compute_LD/0_EUR_UKBB/submitjobs.sh"))




a <- paste0("#!/usr/bin/env bash
#$ -N ukbb_ld_eur_merge
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

#mkdir  /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/mergeplink
#
#for i in {1..22}
#do
#echo /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/EUR/HDL/ref_chr${i} >> /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/mergeplink/EUR_HDL.txt
#done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \\
--keep-allele-order \\
--merge-list /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/mergeplink/EUR_HDL.txt \\
--make-bed \\
--threads 3 \\
--out /dcs04/nilanjan/data/jzhang2/MEPRS/revision1/LD_reference/GLGC/ref_geno/EUR/HDL/allchr

")

writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum/0_compute_LD/0_EUR_UKBB/ukbb_ld_eur_merge.sh"))

cd /dcs04/nilanjan/data/jzhang2/MEPRS/codes/revision1/LD_reference/GLGC/multi-lassosum/0_compute_LD/0_EUR_UKBB
qsub -l h=${skipnode} -hold_jid ukbb_ld_eur ukbb_ld_eur_merge.sh

