

library(readr)
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(readr)
library(cowplot)

My_Theme = theme(
  panel.background = element_blank(), 
  title = element_text(size = 8.5),
  text = element_text(size = 8),
  strip.text.x = element_text(size = 8),
  strip.text.y = element_text(size = 8), 
  axis.title.x = element_text(size = 8.5),
  axis.text.x = element_text(size = 8),
  axis.title.y = element_text(size = 8.5),
  axis.text.y = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
)


######################################################
######################################################
## GLGC


## 1. Regenerate previous results


m2="multi-lassosum_train_from_single_eth"

load("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/GLGC/clean_results/glgc.prediction.result.summary.rdata")
load("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/GLGC/clean_results/weighted_prs_glgc.rdata")
final_result <- final_result[final_result$eth != "EUR",]
prediction.result <- prediction.result[,c(-3,-6,-7,-8)]
for (i in 1:nrow(final_result)) {
  prediction.result$result[(prediction.result$eth == final_result$eth[i]) & (prediction.result$trait == final_result$trait[i]) & (prediction.result$method == "Weighted PRS (CT)")] <- final_result$r2[i]
}
colnames(prediction.result) <- c("eth.vec","t_vec","r2.vec","method_vec")
prediction.result$method_vec <- as.character(prediction.result$method_vec)
prediction.result$method_vec[prediction.result$method_vec == "Best EUR SNP (CT)"] <- "EUR CT"
prediction.result$method_vec[prediction.result$method_vec == "Best EUR PRS (LDpred2)"] <- "EUR LDpred2"
prediction.result$method_vec[prediction.result$method_vec == "Weighted PRS (LDpred2)"] <- "weighted LDpred2"
prediction.result$method_vec[prediction.result$method_vec == "Weighted PRS (CT)"] <- "weighted CT"
prediction.result$method_vec[prediction.result$method_vec == "CT-SLEB (five ancestries)"] <- "CT-SLEB"
prediction.result$method_vec[prediction.result$method_vec == "PRS-CSx (five ancestries)"] <- "PRS-CSx"
prediction.result <- prediction.result[! prediction.result$method_vec %in% c("MEBayes","MELasso","PolyPred+","XPASS",
                                                                             "CT-SLEB (two ancestries)","PRS-CSx (two ancestries)"),]

dat <- read_tsv(paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/GLGC/clean_results/table_",m2,".txt"))
dat <- dat[! dat$method_vec %in% c("weighted lassosum2 (2)",
                                   "multi-ethnic lasso (best)","multi-ethnic lasso (lasso)"),]
dat$method_vec[dat$method_vec == "multi-ethnic lasso (super learning)"] <- "PROSPER"
dat$method_vec[dat$method_vec == "weighted lassosum2 (all ethnic)"] <- "weighted lassosum2"
prediction.result <- rbind(prediction.result,dat)


## 2. Update new bigsnpr results (LDpred2 and lassosum2) 
## 2.1 Remove the previous LDpred2 and lassosum2 results
prediction.result <- prediction.result[! prediction.result$method_vec %in% c("LDpred2","weighted LDpred2","EUR LDpred2","lassosum2", "weighted lassosum2", "EUR lassosum2"),]
## 2.2 Update to the new results
## 2.2.1 lassosum2
tmp <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/GLGC/clean_results/GLGC_lassosum2_sl_rev2.txt")
tmp <- reshape2::melt(tmp[,1:5])
colnames(tmp) <- c("t_vec","eth.vec","method_vec","r2.vec")
tmp$method_vec <- as.character(tmp$method_vec)
tmp$method_vec[tmp$method_vec == "lassosum2_wt"] <- "weighted lassosum2"
tmp$method_vec[tmp$method_vec == "lassosum2_EUR"] <- "EUR lassosum2"
prediction.result <- rbind(prediction.result, tmp)
## 2.2.2 ldpred2
tmp <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/GLGC/clean_results/GLGC_ldpred2_rev2.txt")
tmp <- reshape2::melt(tmp[,1:5])
colnames(tmp) <- c("t_vec","eth.vec","method_vec","r2.vec")
tmp$method_vec <- as.character(tmp$method_vec)
tmp$method_vec[tmp$method_vec == "ldpred2"] <- "LDpred2"
tmp$method_vec[tmp$method_vec == "ldpred2_wt"] <- "weighted LDpred2"
tmp$method_vec[tmp$method_vec == "ldpred2_EUR"] <- "EUR LDpred2"
prediction.result <- rbind(prediction.result, tmp)




library(RColorBrewer)

Single_ancestry_method <- c("CT","LDpred2","lassosum2")
EUR_PRS_based_method <- c("EUR CT","EUR LDpred2","EUR lassosum2")
Multi_ancestry_method_0 <- c("weighted CT","weighted LDpred2","weighted lassosum2")
Multi_ancestry_method_1 <- c("PRS-CSx","CT-SLEB")
Multi_ancestry_method_2 <- c( "PROSPER")

prediction.result$method_vec <- factor(prediction.result$method_vec,
                                       levels = c(Single_ancestry_method, EUR_PRS_based_method, Multi_ancestry_method_0,Multi_ancestry_method_1, Multi_ancestry_method_2))


prediction.result <- prediction.result[! prediction.result$eth.vec %in% c("EUR","AMR"),]
prediction.result <- prediction.result[prediction.result$t_vec != "nonHDL",]
prediction.result <- prediction.result[prediction.result$method_vec != "",]



single.color =  brewer.pal(9, "Blues")[c(3,5,7)]
EUR.color = brewer.pal(9, "Reds")[c(3,5,7)]
multi.color0 = brewer.pal(9, "Greens")[c(4,6,8)] 
multi.color1 = brewer.pal(9, "Oranges")[c(4,6)]
multi.color2 = brewer.pal(9, "Purples")[c(7)] 

colour = c(single.color, EUR.color, multi.color0, multi.color1,multi.color2)

source("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Results/codes_plotting/theme_publication.R")

col_df = tibble(
  colour = colour,
  method_vec = factor(c(Single_ancestry_method, EUR_PRS_based_method, Multi_ancestry_method_0, Multi_ancestry_method_1,Multi_ancestry_method_2),
                      levels = c(Single_ancestry_method, EUR_PRS_based_method, Multi_ancestry_method_0, Multi_ancestry_method_1,Multi_ancestry_method_2)),
  category = case_when(method_vec%in%Single_ancestry_method ~ "Single ancestry method",
                       method_vec%in%EUR_PRS_based_method ~ "EUR PRS based method",
                       method_vec%in%Multi_ancestry_method_0 ~ "Multi ancestry method\n(weighted PRS)",
                       method_vec%in%Multi_ancestry_method_1 ~ "Multi ancestry method\n(existing methods)",
                       method_vec%in%Multi_ancestry_method_2 ~ "PROSPER"
  )
) %>% 
  mutate(category = factor(category,levels = c("Single ancestry method",
                                               "EUR PRS based method",
                                               "Multi ancestry method\n(weighted PRS)",
                                               "Multi ancestry method\n(existing methods)",
                                               "PROSPER")))



prediction.result = prediction.result %>% 
  left_join(col_df)
getLegend <- function(p) {
  g <- ggplotGrob(p)
  k <- which(g$layout$name=="guide-box")
  g$grobs[[k]]
}

run_plot = function(filler, values) {
  values = col_df %>% 
    filter(category %in% filler)
  labels = values %>% 
    pull(method_vec)
  values = values %>% pull(colour)
  names(values) = labels
  ggplot(
    prediction.result %>% 
      filter(category %in% filler),
    aes(x= method_vec, y=r2.vec))+
    geom_bar(aes(fill = method_vec),
             stat="identity",
             position = position_dodge())+
    #geom_point(aes(color=method_vec))+
    theme_Publication()+
    ylab("R2")+
    xlab("Sample Size")+
    labs(fill = filler)+
    scale_fill_manual(values = values)+
    My_Theme
}


legs = lapply(sort(unique(col_df$category)), run_plot)

legs = lapply(legs, getLegend)
p.leg = plot_grid(NULL,
                  legs[[1]],legs[[2]], legs[[3]],
                  legs[[4]],
                  legs[[5]],
                  NULL,
                  align="v",ncol=1,
                  rel_heights=c(0.3,
                                0.7,0.7,0.7,
                                0.6,
                                0.5,
                                0.3))
# p.leg
################################################


prediction.result.sub = prediction.result

p.null <- ggplot(prediction.result.sub,aes(x= method_vec, y=r2.vec))+
  geom_bar(aes(fill=method_vec),
           stat="identity",
           position = position_dodge())+
  theme_Publication()+
  ylab(bquote(''
              *R^2*'                                  '
              *R^2*'                                  '
              *R^2*'                                  '
              *R^2*''))+
  xlab(NULL)+
  labs(fill = "Method")+
  facet_grid(vars(t_vec),vars(eth.vec),scales="free")+
  scale_fill_manual(values = colour) +
  My_Theme+
  scale_x_discrete(labels = NULL, breaks = NULL)+
  theme(legend.position = "none") 


################################################

tmp <- data.frame()
p.num <- ggplot(tmp) + geom_point() + theme_void()


library(ggpubr)

p <- ggarrange(ggarrange(
                         p.num,p.num,p.num,p.num,p.num,
                         ncol = 1,
                         labels = c(NA,"a", "b", "c","d"),
                         heights = c(0.05, 0.95/4, 0.95/4, 0.95/4, 0.95/4)),
               p.null,
               p.leg,
               nrow = 1, 
               labels = c(NA,NA),
               widths = c(0.03,1-0.2-0.03,0.2)
)

ggsave(filename=paste0("GLGC.pdf"),
       plot=p, device="pdf",
       path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision2/updated_figures/", 
       width=180, height=160, units="mm", dpi=320)


################################################

a <- prediction.result[,c(1:3,6,4)]
a$category <- as.character(a$category)
a$category[a$category == "Multi ancestry method\n(weighted PRS)"] <- "Multi ancestry method (weighted PRS)"
a$category[a$category == "Multi ancestry method\n(existing methods)"] <- "Multi ancestry method (existing methods)"

a$category <- factor(a$category, levels = c("PROSPER" ,
                                            "Single ancestry method", "EUR PRS based method",
                                            "Multi ancestry method (weighted PRS)","Multi ancestry method (existing methods)"))
a <- a[order(a$category),]
a <- a[a$eth.vec != "EUR",]
write_csv(a, paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision2/cleaned_tab/GLGC.csv"))


