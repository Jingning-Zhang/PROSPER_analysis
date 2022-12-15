
rm(list=ls())

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/snps")

for (ethnic in c("EUR","AFR","AMR")){
  for (trait in c("bmi","height")){
    system(paste0("awk 'NR>1{print $1}'",
                  " /dcs04/nilanjan/data/jzhang2/summary_data/allofus_cleaned/",ethnic,"/",trait,".txt",
                  " > /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/allofus/snps/",ethnic,"_",trait,"_rsid.txt"))
    print(paste0(ethnic,"_",trait))
  }
}

