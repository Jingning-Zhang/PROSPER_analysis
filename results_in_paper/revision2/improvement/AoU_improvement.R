## AoU
library(readr)
dat <- read_csv(paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision2/cleaned_tab/AoU.csv"))

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

# weighted LDpred2            CT-SLEB            PRS-CSx weighted lassosum2          lassosum2        weighted CT 
# 0.6602678          0.7653043          0.9128994          0.9573457          6.7992354          7.4683932 
# EUR lassosum2        EUR LDpred2            LDpred2             EUR CT                 CT 
# 8.3941859         10.4896749         10.6382289         15.3454797         17.6311570 

mean(res$`PRS-CSx`)
# 0.9128994
mean(res$`CT-SLEB`)
# 0.7653043


library(reshape2)
tab <- melt(res)
tab <- tab[tab$variable != "PROSPER",]
tab$value <- tab$value*100
write_csv(tab, "/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision2/cleaned_tab/AoU_improvement.csv")



