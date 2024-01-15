
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



switchflex <- function(x){
  y <- x
  if(x %in% c("Best EUR SNP (CT)",
              "Weighted PRS (CT)",
              "CT-SLEB (five ancestries)",
              "PRS-CSx (five ancestries)"
              
  )){
    y <- switch (x,
                 "Best EUR SNP (CT)" = "EUR CT",
                 "Weighted PRS (CT)" = "weighted CT",
                 "CT-SLEB (five ancestries)" = "CT-SLEB",
                 "PRS-CSx (five ancestries)" = "PRS-CSx",
                 
    )
  }
  return(y)
}


# GA=1

for(GA in 1:5){
  

##############################
## Results of other methods

load("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Simulation/simu_result_rep.rdata")
result_rep <- result_rep[result_rep$rep_vec %in% 1:3,]
prediction.result <- result_rep %>% 
  group_by(eth_vec,l_vec,ga_vec,m_vec,method_vec) %>% 
  summarise(r2.vec=mean(r2_vad))

prediction.result$method_vec <- as.character(prediction.result$method_vec)

tmp <- prediction.result$method_vec; tmp <- unlist(sapply(tmp, function(x){ switchflex(x) })); names(tmp) <- NULL; prediction.result$method_vec <- tmp
# prediction.result_ct <- prediction.result[prediction.result$method_vec=="CT",]
# View(prediction.result[prediction.result$method_vec=="weighted CT",])
# prediction.result_nonct <- prediction.result[!prediction.result$method_vec=="CT",]
# colnames(prediction.result_nonct) <- c("eth_vec","l_vec","m_vec","ga_vec","method_vec","r2.vec")
# prediction.result <- rbind(prediction.result_ct, prediction.result_nonct)

prediction.result <- prediction.result[prediction.result$m_vec %in% 1:4,]

##############################
## Results of new lassosum2 

tmp <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Simulation/Simulation_lassosum2_sl_rev2_rep_ave.txt")
tmp <- reshape2::melt(tmp, id=c("setting","rho","size","ethnic"))
tmp <- tmp[tmp$setting==GA,]
colnames(tmp) <- c("ga_vec","l_vec","m_vec","eth_vec","method_vec","r2.vec")
tmp$method_vec <- as.character(tmp$method_vec)
tmp$method_vec[tmp$method_vec == "lassosum2_wt"] <- "weighted lassosum2"
tmp$method_vec[tmp$method_vec == "lassosum2_EUR"] <- "EUR lassosum2"

prediction.result <- rbind(data.frame(prediction.result), 
                           data.frame(tmp[tmp$method_vec %in% c("lassosum2","weighted lassosum2","EUR lassosum2","PROSPER"),]))

##############################
## Results of new LDpred2

tmp <- read_tsv("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Simulation/Simulation_ldpred2_rev2_rep_ave.txt")
tmp <- reshape2::melt(tmp[,1:7], id=c("setting","rho","size","ethnic"))
tmp <- tmp[tmp$setting==GA,]
colnames(tmp) <- c("ga_vec","l_vec","m_vec","eth_vec","method_vec","r2.vec")
tmp$method_vec <- as.character(tmp$method_vec)
tmp$method_vec[tmp$method_vec == "ldpred2"] <- "LDpred2"
tmp$method_vec[tmp$method_vec == "ldpred2_wt"] <- "weighted LDpred2"
tmp$method_vec[tmp$method_vec == "ldpred2_EUR"] <- "EUR LDpred2"

prediction.result <- rbind(data.frame(prediction.result),
                           data.frame(tmp[tmp$method_vec %in% c("LDpred2","weighted LDpred2","EUR LDpred2"),]))


##############################
## Plot

Single_ethnic_method <- c("CT","LDpred2","lassosum2")
EUR_PRS_based_method <- c("EUR CT","EUR LDpred2","EUR lassosum2")
Multi_ethnic_method_0 <- c("weighted CT","weighted LDpred2","weighted lassosum2")
Multi_ethnic_method_1 <- c("PRS-CSx","CT-SLEB")
Multi_ethnic_method_2 <- c("PROSPER")

single.color =  brewer.pal(9, "Blues")[c(3,5,7)]
EUR.color = brewer.pal(9, "Reds")[c(3,5,7)]
multi.color0 = brewer.pal(9, "Greens")[c(4,6,8)] ## weighted 
multi.color1 = brewer.pal(9, "Oranges")[c(4,6)] ## exsiting
multi.color2 = brewer.pal(9, "Purples")[c(7)] ## PROSPER

colour = c(single.color, EUR.color, multi.color0, multi.color1, multi.color2)

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
                             levels = c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0, Multi_ethnic_method_1, Multi_ethnic_method_2))) %>% 
  mutate(ga_arc = case_when(ga_vec==1 ~"Strong negative selection (fixed common-SNP heritability)",
                            ga_vec==2 ~"Strong negative selection (fixed per-SNP heritability)",
                            ga_vec==3 ~"Less correlation",
                            ga_vec==4 ~"No negative selection",
                            ga_vec==5 ~"Mild negative selection"
  ))


col_df = tibble(
  colour = colour,
  method_vec = factor(c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0, Multi_ethnic_method_1, Multi_ethnic_method_2),
                      levels = c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0, Multi_ethnic_method_1, Multi_ethnic_method_2)),
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

library(cowplot)


p.null <- list()
for(m in 1:4){
  
  prediction.result.sub <- prediction.result %>% 
    filter(ga_vec==GA &
             eth_vec!="EUR" &
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
    facet_grid(vars(cau_vec),vars(eth_vec))+
    #scale_fill_nejm()+
    scale_fill_manual(values = colour) +
    My_Theme +
    ggtitle(title)+
    theme(legend.position = "none")
  
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

library(ggpubr)

  p <- ggarrange(ggarrange(p.null[[1]],p.null[[3]],
                           nrow = 2, 
                           labels = c("a", "b"),
                           heights = c(0.5,0.5)),
                 p.leg,
                 ncol = 2, widths = c(0.8,0.2))
  ggsave(filename=paste0("GA_",GA,".pdf"),
         plot=p, device="pdf",
         path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/revision2/updated_figures/", 
         width=180, height=180, units="mm", dpi=320)
  
  
  
  a <- prediction.result[prediction.result$ga_vec==GA,c(1,7,8,6,11,5,9)]
  a$category <- as.character(a$category)
  a$category[a$category == "Multi ethnic method\n(weighted PRS)"] <- "Multi ethnic method (weighted PRS)"
  a$category[a$category == "Multi ethnic method\n(existing methods)"] <- "Multi ethnic method (existing methods)" 
  a$category <- factor(a$category, levels = c("PROSPER" , 
                                              "Single ethnic method", "EUR PRS based method",
                                              "Multi ethnic method (weighted PRS)","Multi ethnic method (existing methods)"))
  a$cau_vec <- as.character(a$cau_vec)
  a$cau_vec[a$cau_vec == "p_cau =\n0.01"] <- "p_cau = 0.01"
  a$cau_vec[a$cau_vec == "p_cau =\n0.001"] <- "p_cau = 0.001"
  a$cau_vec[a$cau_vec == "p_cau =\n5E-04"] <- "p_cau = 5E-04" 
  
  a <- a[order(a$eth_vec),]
  a <- a[order(a$sample_size),]
  a <- a[order(a$cau_vec),]
  a <- a[order(a$category),]
  a <- a[a$eth_vec != "EUR",]
  write_csv(a, paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/revision2/cleaned_tab/Simulation_GA_",GA,".csv"))
  
  

}
