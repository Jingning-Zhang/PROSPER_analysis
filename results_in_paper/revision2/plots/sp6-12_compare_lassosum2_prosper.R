
################################################
################################################

## GLGC

################################################
################################################


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


tmp <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/GLGC/clean_results/GLGC_lassosum2_sl_rev2.txt")
tmp <- reshape2::melt(tmp)
colnames(tmp) <- c("t_vec","eth.vec","method_vec","r2.vec")
tmp$method_vec <- as.character(tmp$method_vec)
prediction.result <- tmp


prediction.result <- prediction.result[! prediction.result$eth.vec %in% c("EUR","AMR"),]
prediction.result <- prediction.result[prediction.result$t_vec != "nonHDL",]
prediction.result <- prediction.result[!prediction.result$method_vec %in% c("lassosum2_EUR"),]
prediction.result$method_vec[prediction.result$method_vec == "lassosum2_sl" ] <- "advanced weighted lassosum2"
prediction.result$method_vec[prediction.result$method_vec == "lassosum2_wt"] <- "weighted lassosum2"

# prediction.result <- prediction.result[prediction.result$method_vec != "lassosum2",]

library(RColorBrewer)

Single_ethnic_method <- c("lassosum2", "weighted lassosum2", "advanced weighted lassosum2")
Multi_ethnic_method <- c( "PROSPER" )

prediction.result$method_vec <- factor(prediction.result$method_vec,
                                       levels = c(Single_ethnic_method, Multi_ethnic_method))

single.color =  brewer.pal(9, "Blues")[c(4,6,8)]
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

p

ggsave(filename=paste0("GLGC_sl.pdf"),
       plot=p, device="pdf",
       path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision2/lassosum_sl", 
       width=180, height=160, units="mm", dpi=320)


################################################

prediction.result <- prediction.result[prediction.result$method_vec %in% c("advanced weighted lassosum2","PROSPER"),]
prediction.result <- prediction.result[,c(2,1,4,3)]

prediction.result <- prediction.result[order(prediction.result$eth.vec),]
prediction.result <- prediction.result[order(prediction.result$t_vec),]


prediction.result <- prediction.result[prediction.result$eth.vec != "EUR",]
write_csv(prediction.result, paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision2/GLGC_sl.csv"))




################################################
################################################

## AoU

################################################
################################################




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


tmp <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/AoU/clean_results/AoU_lassosum2_sl_rev2.txt")
tmp <- reshape2::melt(tmp)
colnames(tmp) <- c("t_vec","eth.vec","method_vec","r2.vec")
tmp$method_vec <- as.character(tmp$method_vec)
prediction.result <- tmp


prediction.result <- prediction.result[! prediction.result$eth.vec %in% c("EUR","AMR"),]
prediction.result$t_vec[prediction.result$t_vec == "height"] <- "Height"
prediction.result$t_vec[prediction.result$t_vec == "bmi"] <- "BMI"
prediction.result <- prediction.result[!prediction.result$method_vec %in% c("lassosum2_EUR"),]
prediction.result$method_vec[prediction.result$method_vec == "lassosum2_sl" ] <- "advanced weighted lassosum2"
prediction.result$method_vec[prediction.result$method_vec == "lassosum2_wt"] <- "weighted lassosum2"



# prediction.result <- prediction.result[prediction.result$method_vec != "lassosum2",]

library(RColorBrewer)

Single_ethnic_method <- c("lassosum2", "weighted lassosum2", "advanced weighted lassosum2")
Multi_ethnic_method <- c( "PROSPER" )

prediction.result$method_vec <- factor(prediction.result$method_vec,
                                       levels = c(Single_ethnic_method, Multi_ethnic_method))

single.color =  brewer.pal(9, "Blues")[c(4,6,8)]
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
              *R^2*'                                                         '
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
  heights = c(0.1, 0.9/2, 0.9/2)),
  p.null,
  p.leg,
  nrow = 1, 
  labels = c(NA,NA),
  widths = c(0.03,1-0.3-0.03,0.3)
)

p

ggsave(filename=paste0("AoU_sl.pdf"),
       plot=p, device="pdf",
       path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision2/lassosum_sl/", 
       width=160, height=140, units="mm", dpi=320)


################################################


prediction.result <- prediction.result[prediction.result$method_vec %in% c("advanced weighted lassosum2","PROSPER"),]
prediction.result <- prediction.result[,c(2,1,4,3)]

prediction.result <- prediction.result[order(prediction.result$eth.vec),]
prediction.result <- prediction.result[order(prediction.result$t_vec),]


prediction.result <- prediction.result[prediction.result$eth.vec != "EUR",]
write_csv(prediction.result, paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision2/AoU_sl.csv"))




################################################
################################################
################################################
## simulations
################################################
################################################
################################################



library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(readr)
library(ggpubr)
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


# GA=1


for(GA in 1:5){
  
  tmp <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Simulation/Simulation_lassosum2_sl_rev2_rep_ave.txt")
  tmp <- reshape2::melt(tmp[,c(1:5,7:9)], id=c("setting","rho","size","ethnic"))
  prediction.result <- tmp[tmp$setting==GA,]
  colnames(prediction.result) <- c("ga_vec","l_vec","m_vec","eth.vec","method_vec","r2.vec")
  prediction.result$method_vec <- as.character(prediction.result$method_vec)
  
  prediction.result$method_vec[prediction.result$method_vec == "lassosum2_wt"] <- "weighted lassosum2"
  prediction.result$method_vec[prediction.result$method_vec == "lassosum2_sl" ] <- "advanced weighted lassosum2"
  
  
  Single_ethnic_method <- c("lassosum2", "weighted lassosum2", "advanced weighted lassosum2")
  Multi_ethnic_method <- c( "PROSPER" )
  
  
  single.color =  brewer.pal(9, "Blues")[c(4,6,8)]
  multi.color = brewer.pal(9, "Purples")[c(7)] 
  
  colour = c(single.color, multi.color)
  
  source("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Results/codes_plotting/theme_publication.R")
  
  
  prediction.result = prediction.result %>% 
    mutate(cau_vec = case_when(
      l_vec== 1 ~ paste0("p_cau =\n0.01"),
      l_vec== 2 ~ paste0("p_cau =\n0.001"),
      l_vec== 3 ~ paste0("p_cau =\n5E-04"),
      l_vec== 4 ~ paste0("p_cau =\n0.05"),
      l_vec== 5 ~ paste0("p_cau =\n0.1")
    ),
    sample_size = case_when(
      m_vec == 1~ "15000",
      m_vec == 2~ "45000",
      m_vec == 3~ "80000",
      m_vec == 4~ "100000"
    )) %>% 
    mutate(cau_vec = factor(cau_vec,
                            levels = c("p_cau =\n0.1",
                                       "p_cau =\n0.05",
                                       "p_cau =\n0.01",
                                       "p_cau =\n0.001",
                                       "p_cau =\n5E-04")),
           sample_size = factor(sample_size,
                                levels = c("15000","45000","80000","100000")),
           method_vec = factor(method_vec,
                               levels = c(Single_ethnic_method, Multi_ethnic_method))) %>% 
    mutate(ga_arc = case_when(ga_vec==1 ~"Strong negative selection (fixed common-SNP heritability)",
                              ga_vec==2 ~"Strong negative selection (fixed per-SNP heritability)",
                              ga_vec==3 ~"Less correlation",
                              ga_vec==4 ~"No negative selection",
                              ga_vec==5 ~"Mild negative selection"
    ))
  
  
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
      prediction.result.sub %>% 
        filter(category %in% filler),
      aes(x= sample_size,y=r2.vec,
          group=method_vec))+
      geom_bar(aes(fill = method_vec),
               stat="identity",
               position = position_dodge())+
      #geom_point(aes(color=method_vec))+
      theme_Publication()+
      My_Theme+
      xlab("Sample Size")+
      labs(fill = filler)+
      scale_fill_manual(values = values)
  }
  
  
  p.null <- list()
  for(m in 1:4){
    
    prediction.result.sub <- prediction.result %>% 
      filter(ga_vec==GA &
               eth.vec!="EUR" &
               m_vec ==m)
    
    title = prediction.result.sub$ga_arc[1]
    
    p.null[[m]] <- ggplot(prediction.result.sub,aes(x= sample_size,y=r2.vec,group=method_vec))+
      geom_bar(aes(fill=method_vec),
               stat="identity",
               position = position_dodge())+
      #geom_point(aes(color=method_vec))+
      theme_Publication()+
      ylab(bquote(''
                  *R^2*'                  '
                  *R^2*'                  '
                  *R^2*''))+
      xlab("Sample Size")+
      labs(fill = "Method")+
      facet_grid(vars(cau_vec),vars(eth.vec))+
      #scale_fill_nejm()+
      scale_fill_manual(values = colour) +
      My_Theme +
      ggtitle(title)+
      theme(legend.position = "none")
    
  }
  
  
  
  legs = lapply(sort(unique(col_df$category)), run_plot)
  
  legs = lapply(legs, getLegend)
  p.leg = plot_grid(legs[[1]],legs[[2]],
                    align="v",ncol=1,
                    rel_heights=c(0.3,
                                  0.3),
                    rel_widths=c(0.1,0.1))
  
  
  p <- ggarrange(ggarrange(p.null[[1]],p.null[[3]],
                           nrow = 2,
                           labels = c("a", "b"),
                           heights = c(0.5,0.5)),
                 p.leg,
                 ncol = 2, widths = c(0.72,0.28))
  ggsave(filename=paste0("simulation_GA_",GA,".pdf"),
         plot=p, device="pdf",
         path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision2/lassosum_sl/",
         width=180, height=200, units="mm", dpi=320)
  
}



##################################################
## clean source data

tmp <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Simulation/Simulation_lassosum2_sl_rev2_rep_ave.txt")
tmp <- reshape2::melt(tmp[,c(1:5,7:9)], id=c("setting","rho","size","ethnic"))
prediction.result <- tmp
colnames(prediction.result) <- c("ga_vec","l_vec","m_vec","eth.vec","method_vec","r2.vec")
prediction.result$method_vec <- as.character(prediction.result$method_vec)

prediction.result$method_vec[prediction.result$method_vec == "lassosum2_wt"] <- "weighted lassosum2"
prediction.result$method_vec[prediction.result$method_vec == "lassosum2_sl" ] <- "advanced weighted lassosum2"

Single_ethnic_method <- c("lassosum2", "weighted lassosum2", "advanced weighted lassosum2")
Multi_ethnic_method <- c( "PROSPER" )


prediction.result = prediction.result %>% 
  mutate(cau_vec = case_when(
    l_vec== 1 ~ paste0("p_cau = 0.01"),
    l_vec== 2 ~ paste0("p_cau = 0.001"),
    l_vec== 3 ~ paste0("p_cau = 5E-04"),
    l_vec== 4 ~ paste0("p_cau = 0.05"),
    l_vec== 5 ~ paste0("p_cau = 0.1")
  ),
  sample_size = case_when(
    m_vec == 1~ "15000",
    m_vec == 2~ "45000",
    m_vec == 3~ "80000",
    m_vec == 4~ "100000"
  ),
  ga_arc = case_when(ga_vec==1 ~ "Strong negative selection (fixed common-SNP heritability)",
                              ga_vec==2 ~"Strong negative selection (fixed per-SNP heritability)",
                              ga_vec==3 ~"Less correlation",
                              ga_vec==4 ~"No negative selection",
                              ga_vec==5 ~"Mild negative selection"
  ) ) %>% 
  mutate(cau_vec = factor(cau_vec,
                          levels = c("p_cau = 0.1",
                                     "p_cau = 0.05",
                                     "p_cau = 0.01",
                                     "p_cau = 0.001",
                                     "p_cau = 5E-04")),
         sample_size = factor(sample_size,
                              levels = c("15000","45000","80000","100000")),
         method_vec = factor(method_vec,
                             levels = c(Single_ethnic_method, Multi_ethnic_method)),
         ga_arc = factor(ga_arc, 
                         levels = c("Strong negative selection (fixed common-SNP heritability)",
                                    "Strong negative selection (fixed per-SNP heritability)",
                                    "Less correlation",
                                    "No negative selection",
                                    "Mild negative selection"))
         ) 

prediction.result <- prediction.result[prediction.result$method_vec %in% c("advanced weighted lassosum2","PROSPER"),]
prediction.result <- prediction.result[,c(4,7,8,6,5,9)]

prediction.result <- prediction.result[order(prediction.result$eth.vec),]
prediction.result <- prediction.result[order(prediction.result$sample_size),]
prediction.result <- prediction.result[order(prediction.result$cau_vec),]
prediction.result <- prediction.result[order(prediction.result$ga_arc),]

write_csv(prediction.result, 
          paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision2/Simulation_sl.csv"))




