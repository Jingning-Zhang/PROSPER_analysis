## GLGC
library(readr)
dat <- read_csv(paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/GLGC.csv"))

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
# View(res)

res_max_impr <- apply(res[ ,-1:-3], 2, mean)
sort(res_max_impr)

# weighted lassosum2        weighted CT            PRS-CSx            CT-SLEB          lassosum2 
# 0.06343826         0.28650903         0.34164673         0.37722070         0.49818196 
# EUR lassosum2   weighted LDpred2                 CT            LDpred2             EUR CT 
# 0.79001004         0.88895268         0.94907112         1.32172089         2.85228790 
# EUR LDpred2 
# 2.89937058 

mean(res[ (res$ethnicity_vec %in% c("AFR")),][["weighted lassosum2"]])
# 0.1311333
mean(res[ (res$ethnicity_vec %in% c("SAS")),][["weighted lassosum2"]])
# 0.1225099 
mean(res$`PRS-CSx`)
# 0.3416467
mean(res$`CT-SLEB`)
# 0.3772207

library(reshape2)
tab <- melt(res)
tab <- tab[tab$variable != "PROSPER",]
tab$value <- tab$value*100
write_csv(tab, "/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/GLGC_improvement.csv")


