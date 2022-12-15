
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/lassosum2")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/lassosum2/0_3_mergeplink")
for (ethnic in c("EUR","AFR","AMR")){
  for (trait in c("bmi","height")){

    a <- paste0("#!/usr/bin/env bash
#$ -N ",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G,h_fsize=5000G
#$ -m e
#$ -M jzhan218@jhu.edu

for i in {1..22}
do
echo /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/ref_geno/",ethnic,"/",trait,"/ref_chr${i} >> /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/lassosum2/0_3_mergeplink/",ethnic,"_",trait,".txt
done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \\
--keep-allele-order \\
--merge-list /dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/lassosum2/0_3_mergeplink/",ethnic,"_",trait,".txt \\
--make-bed \\
--threads 1 \\
--out /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/ref_geno/",ethnic,"/",trait,"/allchr

")

    writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/allofus/lassosum2/0_3_mergeplink/",ethnic,"_",trait,".sh"))
  }

}


