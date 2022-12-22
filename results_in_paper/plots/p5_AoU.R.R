

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
  strip.text.y = element_text(size = 9), 
  axis.title.x = element_text(size = 8.5),
  axis.text.x = element_text(size = 8),
  axis.title.y = element_text(size = 8.5),
  axis.text.y = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8)
)

######################################################
######################################################
## AoU


m2="multi-lassosum_train_from_single_eth"

load("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/AoU/clean_results/aou.prediction.result.summary.rdata")
load("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/AoU/clean_results/weighted_prs_aou.rdata")
final_result <- final_result[final_result$eth != "EUR",]
final_result$trait[final_result$trait == "height"] <- "Height"
final_result$trait[final_result$trait == "bmi"] <- "BMI"
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
prediction.result$method_vec[prediction.result$method_vec == "CT-SLEB (three ancestries)"] <- "CT-SLEB"
prediction.result$method_vec[prediction.result$method_vec == "PRS-CSx (three ancestries)"] <- "PRS-CSx"
prediction.result <- prediction.result[! prediction.result$method_vec %in% c("MEBayes","MELasso","PolyPred+","XPASS",
                                                                             "CT-SLEB (two ancestries)","PRS-CSx (two ancestries)"),]

dat <- read_tsv(paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/AoU/clean_results/table_",m2,".txt"))
dat <- dat[! dat$method_vec %in% c("weighted lassosum2 (2)",
                                   "multi-ethnic lasso (best)","multi-ethnic lasso (lasso)"),]
dat$method_vec[dat$method_vec == "multi-ethnic lasso (super learning)"] <- "PROSPER"
dat$method_vec[dat$method_vec == "weighted lassosum2 (all ancestry)"] <- "weighted lassosum2"
prediction.result <- rbind(prediction.result,dat)

prediction.result$t_vec[prediction.result$t_vec == "bmi"] <- "BMI"
prediction.result$t_vec[prediction.result$t_vec == "height"] <- "Height"


library(RColorBrewer)

Single_ethnic_method <- c("CT","LDpred2","lassosum2")
EUR_PRS_based_method <- c("EUR CT","EUR LDpred2","EUR lassosum2")
Multi_ethnic_method_0 <- c("weighted CT","weighted LDpred2","weighted lassosum2")
Multi_ethnic_method_1 <- c("PRS-CSx","CT-SLEB")
Multi_ethnic_method_2 <- c( "PROSPER")

prediction.result$method_vec <- factor(prediction.result$method_vec,
                                       levels = c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0,Multi_ethnic_method_1, Multi_ethnic_method_2))


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
  method_vec = factor(c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0, Multi_ethnic_method_1,Multi_ethnic_method_2),
                      levels = c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0, Multi_ethnic_method_1,Multi_ethnic_method_2)),
  category = case_when(method_vec%in%Single_ethnic_method ~ "Single ethnic method",
                       method_vec%in%EUR_PRS_based_method ~ "EUR PRS based method",
                       method_vec%in%Multi_ethnic_method_0 ~ "Multi ethnic method\n(weighted PRS)",
                       method_vec%in%Multi_ethnic_method_1 ~ "Multi ethnic method\n(existing methods)",
                       method_vec%in%Multi_ethnic_method_2 ~ "PROSPER"
  )
) %>% 
  mutate(category = factor(category,levels = c("Single ethnic method",
                                               "EUR PRS based method",
                                               "Multi ethnic method\n(weighted PRS)",
                                               "Multi ethnic method\n(existing methods)",
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
              *R^2*'                                                                     '
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
  p.num,p.num,p.num,
  ncol = 1,
  labels = c(NA,"a", "b"),
  heights = c(0.1, 0.93/2, 0.93/2)),
  p.null,
  p.leg,
  nrow = 1, 
  labels = c(NA,NA),
  widths = c(0.03,1-0.22-0.03,0.22)
)

ggsave(filename=paste0("AoU.pdf"),
       plot=p, device="pdf",
       path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/", 
       width=180, height=140, units="mm", dpi=320)



################################################

a <- prediction.result[,c(1:3,6,4)]
a$category <- as.character(a$category)
a$category[a$category == "Multi ethnic method\n(weighted PRS)"] <- "Multi ethnic method (weighted PRS)"
a$category[a$category == "Multi ethnic method\n(existing methods)"] <- "Multi ethnic method (existing methods)" 

a$category <- factor(a$category, levels = c("PROSPER" , 
                                            "Single ethnic method", "EUR PRS based method",
                                            "Multi ethnic method (weighted PRS)","Multi ethnic method (existing methods)"))
a <- a[order(a$category),]
a <- a[a$eth.vec != "EUR",]
write_csv(a, paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/AoU.csv"))



