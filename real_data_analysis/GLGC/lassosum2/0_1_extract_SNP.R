
rm(list=ls())

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC")
dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/snps")

for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){

    system(paste0("awk 'NR>1{print $1}'",
                  " /dcs04/nilanjan/data/jzhang2/summary_data/GLGC_cleaned/",ethnic,"/",trait,".txt",
                  " > /dcs04/nilanjan/data/jzhang2/MEPRS/real_data_analysis/GLGC/snps/",ethnic,"_",trait,"_rsid.txt"))
    print(paste0(ethnic,"_",trait))
  }
}

