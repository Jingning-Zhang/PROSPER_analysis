

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation/lassosum2/0_mergeplink")
for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){

    a <- paste0("#!/usr/bin/env bash
#$ -N ",ethnic,"
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=5000G
#$ -pe local 8
#$ -m e
#$ -M jzhan218@jhu.edu

for i in {1..22}
do
echo /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/ref_chr${i} >> /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation/lassosum2/0_mergeplink/",ethnic,".txt
done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \\
--keep-allele-order \\
--merge-list /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation/lassosum2/0_mergeplink/",ethnic,".txt \\
--make-bed \\
--threads 8 \\
--out /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/ref_allchr

")

    writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation/lassosum2/0_mergeplink/",ethnic,".sh"))

}



for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){

    a <- paste0("#!/usr/bin/env bash
#$ -N tv_",ethnic,"
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G,h_fsize=5000G
#$ -pe local 8
#$ -m e
#$ -M jzhan218@jhu.edu

for i in {1..22}
do
echo /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/tv_chr${i} >> /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation/lassosum2/0_mergeplink/tv_",ethnic,".txt
done

/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink1/plink \\
--keep-allele-order \\
--merge-list /dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation/lassosum2/0_mergeplink/tv_",ethnic,".txt \\
--make-bed \\
--threads 8 \\
--out /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic,"/tv_allchr

")

    writeLines(a, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Pipeline_simulation/lassosum2/0_mergeplink/tv_",ethnic,".sh"))

}


