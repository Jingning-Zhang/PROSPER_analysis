
library(readr)
library(bigreadr)

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
    print(paste0(trait,"_",ethnic))
    R2 <- as.matrix(read_csv(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_summary/",switchethnic(ethnic),"_",trait,"_lassosum2_validation"), col_types = cols()))
    colnames(R2) <- sort(paste0("s",1:200))
    R2 <- R2[,paste0("s",1:200)]
    best <- which.max(R2)

    df <- fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/Results_to_23andme/",ethnic,"/",trait,"/prs.file"),
                 sep=" ",
                 colClasses = c("character", rep( "numeric", length(R2))))
    df <- df[,c(1,best+1)]
    fwrite2(df, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_best_original/",ethnic,"_",trait,".prs.file"), col.names = F, sep=" ", nThread=1)
  }
}

for (trait in c("any_cvd","depression","heart_metabolic_disease_burden","height","iqb.sing_back_musical_note","migraine_diagnosis","morning_person")){
  for (ethnic in c("AFR","AMR","EAS","EUR","SAS")){
    print(paste0(trait,"_",ethnic))
    df <- fread2(paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_best_original/",ethnic,"_",trait,".prs.file"),
                 sep=" ",
                 colClasses = c("character", "numeric"))
    df <- df[df$V2 != 0,]
    fwrite2(df, paste0("/dcs04/nilanjan/data/23andme/Analysis/JZ/lassosum2/lassosum2_best/",ethnic,"_",trait,".prs.file"), col.names = F, sep=" ", nThread=1)
  }
}


