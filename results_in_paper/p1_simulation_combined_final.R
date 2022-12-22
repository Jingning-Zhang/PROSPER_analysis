
library(ggplot2)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(readr)

for(GA in 1:5){
  


# theme(axis.text = element_text(size = rel(0.8)),
#       axis.title = element_text(size = rel(0.9)),
#       # legend.title = element_text(size = rel(0.8)),
#       # legend.text = element_text(size = rel(0.8)),
#       title = element_text(size = rel(1)),
#       strip.text = element_text(size = rel(0.7))
# )

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

load("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/simulation/prediction.result.summary.rdata")

prediction.result <- prediction.result[,1:6]
prediction.result$method_vec <- as.character(prediction.result$method_vec)

melasso <- read_tsv(paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/simulation/Restuls_GA_",GA,".txt"))

melasso <- melasso[melasso$method_vec %in% c("multi-lassosum_train_from_single_eth (superlearning_combined_R2_cross_ancestry)",
                                             # "multi-lassosum_train_from_single_eth (lasso_combined_R2_cross_ancestry)",
                                             # "multi-lassosum_train_from_single_eth (best_R2_across_ancestry)",
                                             "lassosum2",
                                             "EUR lassosum2",
                                             "Weighted lassosum2"),]


switchflex <- function(x){
  y <- x
  if(x %in% c("Best EUR SNP (CT)",
              "Best EUR PRS (LDpred2)",
              
              "Weighted PRS (CT)",
              "Weighted PRS (LDpred2)",
              "Weighted lassosum2",
              
              "CT-SLEB (five ancestries)",
              
              "PRS-CSx (five ancestries)",
              
              "multi-lassosum_train_from_single_eth (superlearning_combined_R2_cross_ancestry)"
              # "multi-lassosum_train_from_single_eth (lasso_combined_R2_cross_ancestry)"
              # "multi-lassosum_train_from_single_eth (best_R2_across_ancestry)"
  )){
    y <- switch (x,
                 "Best EUR SNP (CT)" = "EUR CT",
                 "Best EUR PRS (LDpred2)"="EUR LDpred2",
                 
                 "Weighted PRS (CT)" = "weighted CT",
                 "Weighted PRS (LDpred2)" = "weighted LDpred2",
                 "Weighted lassosum2" = "weighted lassosum2",
                 
                 "CT-SLEB (five ancestries)" = "CT-SLEB",
                 
                 "PRS-CSx (five ancestries)" = "PRS-CSx",
                 
                 "multi-lassosum_train_from_single_eth (superlearning_combined_R2_cross_ancestry)" = "PROSPER"
                 # "multi-lassosum_train_from_single_eth (lasso_combined_R2_cross_ancestry)" = "PROSPER (lasso)"
                 # "multi-lassosum_train_from_single_eth (best_R2_across_ancestry)" = "PROSPER (best)"
    )
  }
  return(y)
}


Single_ethnic_method <- c("CT","LDpred2","lassosum2")
EUR_PRS_based_method <- c("EUR CT","EUR LDpred2","EUR lassosum2")
Multi_ethnic_method_0 <- c("weighted CT","weighted LDpred2","weighted lassosum2")
Multi_ethnic_method_1 <- c("PRS-CSx","CT-SLEB")
# Multi_ethnic_method_2 <- c("PROSPER (best)","PROSPER (lasso)","PROSPER (super learning)")
# Multi_ethnic_method_2 <- c("PROSPER (lasso)","PROSPER (super learning)")
Multi_ethnic_method_2 <- c("PROSPER")

single.color =  brewer.pal(9, "Blues")[c(3,5,7)]
EUR.color = brewer.pal(9, "Reds")[c(3,5,7)]
multi.color0 = brewer.pal(9, "Greens")[c(4,6,8)] ## weighted 
multi.color1 = brewer.pal(9, "Oranges")[c(4,6)] ## exsiting
# multi.color2 = brewer.pal(9, "Purples")[c(5,7,9)] ## PROSPER
multi.color2 = brewer.pal(9, "Purples")[c(7)] ## PROSPER

colour = c(single.color, EUR.color, multi.color0, multi.color1, multi.color2)

source("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Results/codes_plotting/theme_publication.R")


prediction.result <- rbind(prediction.result,melasso)


tmp <- prediction.result$method_vec; tmp <- unlist(sapply(tmp, function(x){ switchflex(x) })); names(tmp) <- NULL; prediction.result$method_vec <- tmp



prediction.result0 <-  prediction.result



prediction.result <- prediction.result0[prediction.result0$method_vec %in% 
                                          c(Single_ethnic_method, EUR_PRS_based_method, Multi_ethnic_method_0, Multi_ethnic_method_1, Multi_ethnic_method_2),]


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

a <- prediction.result[prediction.result$ga_vec==GA,c(1,7,8,2,6,11,9)]
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
a <- a[order(a$category),]
a <- a[a$eth.vec != "EUR",]
write_csv(a, paste0("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/source_data/simulation/cleaned_tab/Restuls_GA_",GA,".csv"))


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
         path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/simulation/", 
         width=180, height=180, units="mm", dpi=320)
  
  

}
