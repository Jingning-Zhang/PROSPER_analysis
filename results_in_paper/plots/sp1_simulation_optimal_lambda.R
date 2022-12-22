library(MASS)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

set.seed(1)
h2 <- 0.2
m = 1000
ncausal <- floor(0.05*m)
  
N_vec <- seq(from=10,to=500,by=20)
lambda_best <- numeric()
for(i in 1:length(N_vec)){
  print(i)
  N = N_vec[i]
  X = mvrnorm(N, rep(0,m), ar1_cor(m,0.4))
  beta = numeric(m); beta[sample(m, ncausal)] <- rnorm(ncausal); beta <- matrix(beta, ncol=1)
  y_genet <- X %*% beta
  y_err <- rnorm(N, 0, sqrt(var(y_genet) / h2 * (1-h2)))
  Y <- y_err+y_genet
  
  fit <- glmnet::cv.glmnet(x=as.matrix(X), y=Y, family='gaussian', alpha=1)
  lambda_best[i] <- fit$lambda.min
  
}

dat <- data.frame(N_vec,lambda_best)

library(ggplot2)
source("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Analysis/Results/codes_plotting/theme_publication.R")
p = ggplot(dat, aes(x= N_vec, y=lambda_best))+
  geom_point()+
  theme_Publication()+
  ylab("optimal lambda")+
  xlab("sample size (n)")+
  theme(axis.text = element_text(size = rel(0.9)),
        axis.title = element_text(size = rel(1)),
        title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(0.9)))+
  labs(title = "Optimal lambda in lasso")

dir.create("/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure")
ggsave(filename="simulation_optimal_lambda.pdf",
       plot=p, device="pdf",
       path="/Users/jnz_1/Document/JHU/Research/PRS/MEPRS/Manuscript/figure/", 
       width=120, height=120, units="mm", dpi=320)

# Dimention (p) of design matrix X is 1000, and 5% of predictors are randomly selected to be causal. Correlation struction of X is AR1 with rho=0.4.



