rm(list = ls())

#ethnic='EUR'
#trait='HDL'
for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
    a <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/1_prepare_data_lassosum2/",ethnic,"_",trait,".Rout"))
        if(!("> proc.time()" %in% a)){print(paste0(ethnic,"_",trait,".sh")); print(a[length(a)-1])}
  }
}

rm(list = ls())

#ethnic='EUR'
#trait='HDL'
for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
    a <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/2_run_lassosum2/",ethnic,"_",trait,".Rout"))
    if(!("> proc.time()" %in% a)){print(paste0(ethnic,"_",trait,".sh")); print(a[length(a)-1])}
  }
}



#ethnic='EUR'
#trait='HDL'
for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
    a <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/3_lassosum2_clean_score/",ethnic,"_",trait,".Rout"))
    if(!("> proc.time()" %in% a)){print(paste0(ethnic,"_",trait,".sh"))}
  }
}


for (ethnic in c("EUR","AFR","AMR","EAS","SAS")){
  for (trait in c("HDL","LDL","logTG","nonHDL","TC")){
    a <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/real_data_analysis/GLGC/lassosum2/5_ukbb_tuning+validation_best_EUR/",ethnic,"_",trait,".Rout"))
    if(!("> proc.time()" %in% a)){print(paste0(ethnic,"_",trait,".sh"))}
  }
}