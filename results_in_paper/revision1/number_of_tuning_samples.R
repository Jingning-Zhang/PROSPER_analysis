
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(readr)

  

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


setwd("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Results/")

prediction.result <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision1/number_of_tuning_samples.txt")
colnames(prediction.result) <- c("l_vec", "m_vec","eth.vec","Ntuning_vec","r2.vec")
prediction.result$ga_vec <- 1


prediction.result = prediction.result %>% 
  mutate(cau_vec = case_when(
    l_vec== 1 ~ paste0("p_cau =\n0.01"),
    l_vec== 2 ~ paste0("p_cau =\n0.001"),
    l_vec== 3 ~ paste0("p_cau =\n5E-04")
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
                              levels = c("15000","45000","80000","100000")) ,
         Ntuning_vec = factor(Ntuning_vec, 
                              levels = c("5000","3000","1000","500","300","100"))) %>% 
  mutate(ga_arc = case_when(ga_vec==1 ~"Strong negative selection (fixed common-SNP heritability)",
                            ga_vec==2 ~"Strong negative selection (fixed per-SNP heritability)",
                            ga_vec==3 ~"Less correlation",
                            ga_vec==4 ~"No negative selection",
                            ga_vec==5 ~"Mild negative selection"
  ))


p <- ggplot(prediction.result, aes(x= Ntuning_vec, y=r2.vec,
                              color=sample_size,group=interaction(sample_size,cau_vec,eth.vec))) +
  geom_line()+
  geom_point() +
  facet_grid(vars(cau_vec),vars(eth.vec))+
  ylab(bquote(''
              *R^2*'                                 '
              *R^2*'                                 '
              *R^2*''))+
  xlab("Tuning Sample Size")+
  labs(color = "Training\nsample size")+
  scale_color_manual(values = c("#d7191c","#fdae61","#abdda4","#2b83ba")) +
  My_Theme +
  ggtitle(prediction.result$ga_arc[1]) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        panel.grid.major.x = element_blank() )
  



  ggsave(filename=paste0("number_of_tuning_samples.pdf"),
         plot=p, device="pdf",
         path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision1/", 
         width=180, height=120, units="mm", dpi=320)
  
