## AoU
library(readr)
dat <- read_csv(paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/23andme.csv"))

methods <- unique(dat$method_vec)
ethnicity <- unique(dat$eth.vec)
trait <- unique(dat$t_vec)
res0 <- matrix(nrow=length(ethnicity)*length(trait), ncol=length(methods))

for(r in 1:length(methods)){
  pro_impr <- numeric()
  ethnicity_vec <- character()
  trait_vec <- character()
  t=0
  for (i in 1:length(ethnicity)) {
    for(j in 1:length(trait)){
      t=t+1
      tmp <- dat[(dat$eth.vec == ethnicity[i])&(dat$t_vec == trait[j]),]
      a <- tmp$r2.vec[tmp$method_vec == "PROSPER"]
      b <- tmp$r2.vec[tmp$method_vec == methods[r]]
      pro_impr[t] <- (a-b)/b
      ethnicity_vec[t] <- ethnicity[i]
      trait_vec[t] <- trait[j]
    }
  }
  res0[,r] <- pro_impr
}
colnames(res0) <- methods
res <- tibble(ethnicity_vec,trait_vec, as.data.frame(res0))
View(res)


### continuous traits ###

res_max_impr <- apply(res[ (res$trait_vec %in% c("Heart metabolic disease burden","Height")),-1:-3], 2, mean)
sort(res_max_impr)

# PRS-CSx            CT-SLEB weighted lassosum2   weighted LDpred2        weighted CT 
# 0.1254742          0.2571282          0.3865018          0.6810205          1.1933757 
# lassosum2      EUR lassosum2                 CT            LDpred2             EUR CT 
# 1.8364955          4.5855278          4.6473907          4.8936986          8.1155267 
# EUR LDpred2 
# 9.1107566


## second best
apply(res[(res$ethnicity_vec %in% c("African American","Latino")) & (res$trait_vec %in% c("Heart metabolic disease burden","Height")),-1:-3], 1, which.min)
mean(apply(res[(res$ethnicity_vec %in% c("African American","Latino")) & (res$trait_vec %in% c("Heart metabolic disease burden","Height")),-1:-3], 1, min)) 
# 0.2692183


### binary traits ###

res_max_impr <- apply(res[ ! (res$trait_vec %in% c("Heart metabolic disease burden","Height")),-1:-3], 2, mean)
sort(res_max_impr)

# PRS-CSx            CT-SLEB weighted lassosum2   weighted LDpred2      EUR lassosum2 
# 0.01053899         0.01255363         0.01748548         0.02929965         0.03014897 
# EUR LDpred2        weighted CT             EUR CT          lassosum2                 CT 
# 0.04040979         0.04663348         0.05410601         0.11217767         0.11789154 
# LDpred2 
# 0.12761206 



library(reshape2)
tab <- melt(res)
tab <- tab[tab$variable != "PROSPER",]
tab$value <- tab$value*100
write_csv(tab, "/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/23andMe_improvement.csv")



tab1 <- read_csv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/23andMe_improvement.csv")
tab2 <- read_csv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/GLGC_improvement.csv")
tab3 <- read_csv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/AoU_improvement.csv")
im1 <- tab1$value[ (tab1$variable == "PRS-CSx") & (tab1$ethnicity_vec == "African American") & (tab1$trait_vec  %in% c("Heart metabolic disease burden","Height"))]
im2 <- tab2$value[ (tab2$variable == "PRS-CSx") & (tab2$ethnicity_vec == "AFR")]
im3 <- tab3$value[ (tab3$variable == "PRS-CSx") & (tab3$ethnicity_vec == "AFR")]

mean(c(im1,im2,im3)) # 69.2452
mean(im1) # 32.37764
mean(im2) # 76.65662
mean(im3) # 91.28994


im1 <- tab1$value[ (tab1$variable == "CT-SLEB") & (tab1$ethnicity_vec == "African American") & (tab1$trait_vec  %in% c("Heart metabolic disease burden","Height"))]
im2 <- tab2$value[ (tab2$variable == "CT-SLEB") & (tab2$ethnicity_vec == "AFR")]
im3 <- tab3$value[ (tab3$variable == "CT-SLEB") & (tab3$ethnicity_vec == "AFR")]

mean(c(im1,im2,im3)) # 46.00867
mean(im1) # 53.37563
mean(im2) # 27.06431
mean(im3) # 76.53043


