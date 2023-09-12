

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

prediction.result <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/GLGC/clean_results/GLGC_lassosum2_sl.txt")
prediction.result <- prediction.result[prediction.result$method_vec != "PROSPER (best)", ]
prediction.result$method_vec[prediction.result$method_vec == "PROSPER (sl)"] <- "PROSPER"
prediction.result$method_vec[prediction.result$method_vec == "lassosum2 (best)"] <- "lassosum2"
prediction.result$method_vec[prediction.result$method_vec == "lassosum2 (sl)"] <- "lassosum2 with super learning"

prediction.result <- prediction.result[! prediction.result$eth.vec %in% c("EUR","AMR"),]
prediction.result <- prediction.result[prediction.result$t_vec != "nonHDL",]
prediction.result <- prediction.result[prediction.result$method_vec != "",]


library(RColorBrewer)

Single_ethnic_method <- c("lassosum2", "lassosum2 with super learning")
Multi_ethnic_method <- c( "PROSPER" )

prediction.result$method_vec <- factor(prediction.result$method_vec,
                                       levels = c(Single_ethnic_method, Multi_ethnic_method))

single.color =  brewer.pal(9, "Blues")[c(4,6)]
multi.color = brewer.pal(9, "Purples")[c(7)] 

colour = c(single.color, multi.color)

source("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Results/codes_plotting/theme_publication.R")

col_df = tibble(
  colour = colour,
  method_vec = factor(c(Single_ethnic_method, Multi_ethnic_method),
                      levels = c(Single_ethnic_method, Multi_ethnic_method)),
  category = case_when(method_vec%in%Single_ethnic_method ~ "lassosum2",
                       method_vec%in%Multi_ethnic_method ~ "PROSPER"
  )
) %>% 
  mutate(category = factor(category,levels = c("lassosum2",
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
                  legs[[1]],legs[[2]],
                  NULL,
                  align="v",ncol=1,
                  rel_heights=c(0.3,
                                # 0.7,0.7,0.7,
                                # 0.6,
                                # 0.5,
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
               widths = c(0.03,1-0.3-0.03,0.3)
)

ggsave(filename=paste0("GLGC.pdf"),
       plot=p, device="pdf",
       path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision1/", 
       width=180, height=160, units="mm", dpi=320)


################################################

a <- prediction.result[,c(1:3,6,4)]
a$category <- as.character(a$category)
a$method_vec <- as.character(a$method_vec)

a <- a[order(a$t_vec),]
a <- a[a$eth.vec != "EUR",]
write_csv(a, paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision1/GLGC_sl.csv"))



