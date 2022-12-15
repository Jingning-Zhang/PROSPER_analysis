
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/0_3_mergeplink")

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/mergeplink")
for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){

    a <- paste0("#!/usr/bin/env bash
#$ -N ",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=5000G
#$ -pe local 3
#$ -m e
#$ -M jzhan218@jhu.edu

for i in {1..22}
do
echo /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/",ethnic,"/",trait,"/ref_chr${i} >> /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/mergeplink/",ethnic,"_",trait,".txt
done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \\
--keep-allele-order \\
--merge-list /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/mergeplink/",ethnic,"_",trait,".txt \\
--make-bed \\
--threads 3 \\
--out /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/ref_geno/",ethnic,"/",trait,"/allchr

")
    writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/0_3_mergeplink/",ethnic,"_",trait,".sh"))
  }

}


