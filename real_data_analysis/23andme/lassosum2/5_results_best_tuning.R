
library(readr)

switchethnic <- function(x){
  switch (x,
          "AFR" = "african_american",
          "AMR" = "latino",
          "EAS" = "east_asian",
          "EUR" = "european",
          "SAS" = "south_asian"
  )
}


for (trait in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){
  for (ethnic in c("AFR","AMR","EAS","EUR","SAS")){
    R2 <- as.matrix(read_csv(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_summary/",switchethnic(ethnic),"_",trait,"_lassosum2_validation"), col_types = cols()))
    colnames(R2) <- sort(paste0("s",1:200))
    R2 <- R2[,paste0("s",1:200)]

    best <- which.max(R2)

    para <- read_tsv(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results/para_",ethnic,"_",trait,".txt"),col_types = cols())
    para$NUM <- 1:nrow(para)
    best_para <- para[best,]
    write_tsv(best_para, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_summary/best_para/",ethnic,"_",trait,".txt"))
    print(paste0(ethnic,"-",trait,": delta=",round(best_para$delta,4), "; R2=",round(max(R2),4), "; sparsity=",round(best_para$sparsity,2) ))
  }
}

#[1] "AFR-any_cvd: delta=15.6977; R2=0.5526; sparsity=0.74"
#[1] "AMR-any_cvd: delta=25.1902; R2=0.5922; sparsity=0.86"
#[1] "EAS-any_cvd: delta=100; R2=0.5619; sparsity=0.31"
#[1] "EUR-any_cvd: delta=8.9514; R2=0.6565; sparsity=0.55"
#[1] "SAS-any_cvd: delta=54.2894; R2=0.5348; sparsity=0.46"
#[1] "AFR-depression: delta=37.8978; R2=0.5374; sparsity=0.09"
#[1] "AMR-depression: delta=74.8338; R2=0.569; sparsity=0.53"
#[1] "EAS-depression: delta=37.8978; R2=0.5664; sparsity=0.14"
#[1] "EUR-depression: delta=8.9514; R2=0.6104; sparsity=0.65"
#[1] "SAS-depression: delta=0.5; R2=0.5344; sparsity=0.15"
#[1] "AFR-heart_metabolic_disease_burden: delta=25.1902; R2=0.0134; sparsity=0.25"
#[1] "AMR-heart_metabolic_disease_burden: delta=37.8978; R2=0.0352; sparsity=0.89"
#[1] "EAS-heart_metabolic_disease_burden: delta=100; R2=0.009; sparsity=0.69"
#[1] "EUR-heart_metabolic_disease_burden: delta=8.9514; R2=0.0959; sparsity=0.47"
#[1] "SAS-heart_metabolic_disease_burden: delta=100; R2=0.0052; sparsity=0.74"
#[1] "AFR-height: delta=4.4822; R2=0.0713; sparsity=0.96"
#[1] "AMR-height: delta=15.6977; R2=0.1741; sparsity=0.81"
#[1] "EAS-height: delta=25.1902; R2=0.1233; sparsity=0.92"
#[1] "EUR-height: delta=8.9514; R2=0.2886; sparsity=0.43"
#[1] "SAS-height: delta=15.6977; R2=0.0833; sparsity=0.95"
#[1] "AFR-iqb.sing_back_musical_note: delta=74.8338; R2=0.5449; sparsity=0.59"
#[1] "AMR-iqb.sing_back_musical_note: delta=54.2894; R2=0.5844; sparsity=0.63"
#[1] "EAS-iqb.sing_back_musical_note: delta=100; R2=0.5505; sparsity=0.3"
#[1] "EUR-iqb.sing_back_musical_note: delta=8.9514; R2=0.6712; sparsity=0.59"
#[1] "SAS-iqb.sing_back_musical_note: delta=0.5; R2=0.5671; sparsity=1"
#[1] "AFR-migraine_diagnosis: delta=15.6977; R2=0.5365; sparsity=0.21"
#[1] "AMR-migraine_diagnosis: delta=37.8978; R2=0.5811; sparsity=0.69"
#[1] "EAS-migraine_diagnosis: delta=37.8978; R2=0.5505; sparsity=0.8"
#[1] "EUR-migraine_diagnosis: delta=8.9514; R2=0.6348; sparsity=0.57"
#[1] "SAS-migraine_diagnosis: delta=15.6977; R2=0.5231; sparsity=0.12"
#[1] "AFR-morning_person: delta=8.9514; R2=0.5673; sparsity=0.78"
#[1] "AMR-morning_person: delta=8.9514; R2=0.6195; sparsity=0.91"
#[1] "EAS-morning_person: delta=37.8978; R2=0.5805; sparsity=0.57"
#[1] "EUR-morning_person: delta=8.9514; R2=0.6797; sparsity=0.49"
#[1] "SAS-morning_person: delta=100; R2=0.5498; sparsity=0.68"
