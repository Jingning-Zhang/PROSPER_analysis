
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

source("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Results/codes_plotting/theme_publication.R")
prediction.result <- readRDS("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/23andme/*results_plot/res_Dec.rds") # from /Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/23andme/*results_plot/plot_Dec.R

Single_ethnic_method <- c("CT","lassosum2")
EUR_PRS_based_method <- c("EUR CT","EUR lassosum2")
Multi_ethnic_method_0 <- c("weighted CT","weighted lassosum2")
Multi_ethnic_method_1 <- c("PRS-CSx","CT-SLEB")
Multi_ethnic_method_2 <- c("PROSPER")

prediction.result <- prediction.result[prediction.result$method_vec %in% c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0,Multi_ethnic_method_1, Multi_ethnic_method_2), ]

prediction.result$method_vec <- factor(prediction.result$method_vec,
                                       levels = c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0,Multi_ethnic_method_1, Multi_ethnic_method_2))


colnames(prediction.result) <- c("r2.vec","eth.vec","t_vec","method_vec")

prediction.result <- prediction.result[! prediction.result$eth.vec %in% c("European"),]


single.color =  brewer.pal(9, "Blues")[c(3,7)]
EUR.color = brewer.pal(9, "Reds")[c(3,7)]
multi.color0 = brewer.pal(9, "Greens")[c(4,8)] 
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


prediction.result.sub = prediction.result[!prediction.result$t_vec %in% c("Heart metabolic disease burden","Height"),]

p.null_binary <- ggplot(prediction.result.sub,aes(x= method_vec, y=r2.vec))+
  geom_bar(aes(fill=method_vec),
           stat="identity",
           position = position_dodge())+
  theme_Publication()+
  ylab(bquote('   '
              *AUC*'                             '
              *AUC*'                             '
              *AUC*'                             '
              *AUC*'                             '
              *AUC*''))+
  xlab(NULL)+
  labs(fill = "Method")+
  facet_grid(vars(t_vec),vars(eth.vec),scales="free")+
  scale_fill_manual(values = colour) +
  My_Theme+
  scale_x_discrete(labels = NULL, breaks = NULL)+
  coord_cartesian(ylim = c(0.5, NA))+
  theme(legend.position = "none")


prediction.result.sub = prediction.result[prediction.result$t_vec %in% c("Heart metabolic disease burden","Height"),]

p.null_continuous <- ggplot(prediction.result.sub,aes(x= method_vec, y=r2.vec))+
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


tmp <- data.frame()
p.num <- ggplot(tmp) + geom_point() + theme_void()


library(ggpubr)
p_binary <- ggarrange(ggarrange(
  p.num,p.num,p.num,p.num,p.num,p.num,
  ncol = 1,
  labels = c(NA,"a", "b", "c","d","e"),
  heights = c(0.05, 0.95/5, 0.95/5, 0.95/5, 0.95/5,1/5)),
  p.null_binary,
  p.leg,
  nrow = 1, 
  labels = c(NA,NA),
  widths = c(0.03,1-0.2-0.03,0.2)
)

ggsave(filename=paste0("23andme_binary.pdf"),
       plot=p_binary, device="pdf",
       path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision2/updated_figures/", 
       width=180, height=180, units="mm", dpi=320)


p_continuous <- ggarrange(ggarrange(
  p.num,p.num,p.num,
  ncol = 1,
  labels = c(NA,"a", "b"),
  heights = c(0.05, 0.95/2, 0.95/2)),
  p.null_continuous,
  p.leg,
  nrow = 1, 
  labels = c(NA,NA),
  widths = c(0.03,1-0.22-0.03,0.22)
)


ggsave(filename=paste0("23andme_continuous.pdf"),
       plot=p_continuous, device="pdf",
       path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision2/updated_figures/", 
       width=180, height=140, units="mm", dpi=320)






a <- prediction.result[,c(2,3,1,6,4)]
a$category <- as.character(a$category)
a$category[a$category == "Multi ethnic method\n(weighted PRS)"] <- "Multi ethnic method (weighted PRS)"
a$category[a$category == "Multi ethnic method\n(existing methods)"] <- "Multi ethnic method (existing methods)" 

a$category <- factor(a$category, levels = c("PROSPER" , 
                                            "Single ethnic method", "EUR PRS based method",
                                            "Multi ethnic method (weighted PRS)","Multi ethnic method (existing methods)"))

a <- a[order(a$method_vec),]
a <- a[order(a$category),]
a <- a[a$eth.vec != "EUR",]
write_csv(a, paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision2/cleaned_tab/23andme.csv"))



