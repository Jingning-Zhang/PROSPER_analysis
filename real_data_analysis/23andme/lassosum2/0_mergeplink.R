
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/0_lassosum2/0_mergeplink")
for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){

    a <- paste0("#!/usr/bin/env bash
#$ -N ",ethnic,"_",trait,"
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=5000G
#$ -pe local 8
#$ -m e
#$ -M jzhan218@jhu.edu

for i in {1..22}
do
echo /dcs04/nilanjan/data/23andme/Analysis/JZ/ref_geno/",ethnic,"/",trait,"/ref_chr${i} >> /dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/0_lassosum2/0_mergeplink/",ethnic,"_",trait,".txt
done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \\
--keep-allele-order \\
--merge-list /dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/0_lassosum2/0_mergeplink/",ethnic,"_",trait,".txt \\
--make-bed \\
--threads 8 \\
--out /dcs04/nilanjan/data/23andme/Analysis/JZ/ref_geno/",ethnic,"/",trait,"/allchr

")

    writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/23andme/0_lassosum2/0_mergeplink/",ethnic,"_",trait,".sh"))
  }

}


