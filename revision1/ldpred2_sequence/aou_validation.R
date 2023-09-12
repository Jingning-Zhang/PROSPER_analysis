
# --------------------- Validation ---------------------
rm(list=ls())
library(readr)
library(data.table)
library(mvtnorm)
library(devtools)
library(lavaan)
library(gdata)
library(xtable)
library(MASS) # for the ginv
library(data.table)
library(corpcor) #for pseudoinverse
library(parallel)
library(MendelianRandomization) # for mr_ivw
library(dplyr)
library(R.utils) # for gzip
library(stringr) # for str_detect
library(genio) # a package to facilitate reading and writing genetics data. The focus of this vignette is processing plink BED/BIM/FAM files.
library(data.table)
library(pROC)
library(bigsnpr)
#library(rms)
library(DescTools)
traits = c('height','bmi')
races = c('EUR','AFR','AMR'); K = length(races)

Out = matrix(NA,length(traits),3*length(races))
colnames(Out) = sapply(races,function(x){paste(c('R2','SD','P-value'),x)})
rownames(Out) = traits
n.split = 1

h2_seq <- c(0.3, 0.7, 1, 1.4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
sets <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

pars = matrix(NA,length(traits)*length(races),ncol(sets))
colnames(pars) = colnames(sets)
te = expand.grid(trait=traits,race=races)
rownames(pars) = sapply(1:nrow(te), function(x){paste(te[x,1], te[x,2])})
pars = as.data.frame(pars)




Indx = data.frame(matrix(NA,length(traits),length(races)))
rownames(Indx) = traits; colnames(Indx) = races
R2 = matrix(NA, length(traits), K); colnames(R2) = races; rownames(R2) = traits
for (trait in traits){
  out = matrix(NA,n.split,3*length(races))
  colnames(out) = sapply(races,function(x){paste(c('R2 Adjusted','Regression Coef','P-value'),x)})
  for (race in races){
    writescore_path = paste0('/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/',trait,'/',race,'/ldpred2/')
    validatetable1 = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/',trait,'/tuning+validation/',race,'_tuning.txt'),header=T)
    ids1 = as.character(validatetable1$f.eid)
    validatetable2 = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/',trait,'/tuning+validation/',race,'_validation.txt'),header=T)
    ids2 = as.character(validatetable2$f.eid)
    validatetable = rbind(validatetable1[,1:2], validatetable2[,1:2])
    colnames(validatetable) = c('id','y')
    covariates1 = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/',race,'_tuning.txt'))
    covariates2 = bigreadr::fread2(paste0('/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/covariates/tuning+validation/',race,'_validation.txt'))
    covariates = rbind(covariates1, covariates2)
    colnames(covariates) = c('id','gender','age',paste0('pc',1:10))
    validatetable = merge(validatetable, covariates, by = 'id')
    validatetable = validatetable[complete.cases(validatetable$y),]
    validatetable$id = as.character(validatetable$id)

    for(chr in 1:22){
      temfile = paste0(writescore_path,'score/ldpred2-chr',chr,'.sscore')
      preds = bigreadr::fread2(temfile)
      if(chr==1){
        preds0 = as.matrix(preds[,paste0('SCORE',1:nrow(sets),'_SUM')])
      }else{
        preds0 = as.matrix(preds[,paste0('SCORE',1:nrow(sets),'_SUM')]) + preds0
      }
    }
    preds = cbind(preds[,"IID",drop=F], preds0)
    colnames(preds) = c('id',paste0('prs',1:nrow(sets)))
    preds$id = as.character(preds$id)
    rownames(preds) = preds$id
    
    validatetable = validatetable[validatetable$id %in% rownames(preds),]
    
    valdat = validatetable; rm(validatetable)
    rownames(valdat) = as.character(valdat$id)
    #---------------------------------------#---------------------------------------
    validatetable = valdat[ids1,]
    
    output = matrix(0,nrow(sets),3)
    colnames(output) = c('R2 Adjusted','Regression Coef','P-value')
    rownames(output) = sapply(1:nrow(sets),function(x){paste0('p=',sets[x,1],' h2=', sets[x,2], 'sp=',sets[x,3])})
    
    for(i in 1:nrow(sets)){
      tem = data.frame(id=rownames(preds), prs = preds[,paste0('prs',i)])
      prstable = merge(tem,validatetable,by='id')
      prstable = prstable[complete.cases(prstable),]
      if ((nrow(prstable)>0)&(sum(prstable$prs) != 0)){
        # get residual:
        formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
        fit = lm(formula.prs, data=prstable)
        prstable$res = fit$residuals
        fit = lm(res~prs, data=prstable)
        output[i,'Regression Coef'] = coefficients(fit)['prs']
        output[i,'R2 Adjusted'] = summary(fit)$r.squared #(coefficients(fit)['prs'])^2/var(na.omit(prstable[,'res']))
        output[i,'P-value'] = summary(fit)$coefficients['prs','Pr(>|t|)']
        #print(i)
      }
    }
    Indx[trait,race] = which.max(output[,'R2 Adjusted'])
    tem = data.frame(id=rownames(preds), prs = preds[,paste0('prs',Indx[trait,race])])
    prstable = merge(tem,validatetable,by='id')
    prstable = prstable[complete.cases(prstable),]
    formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
    fit = lm(formula.prs, data=prstable)
    prstable$res = fit$residuals
    fit.train = lm(res~prs, data=prstable)
    print(output[Indx[trait,race],])
    
    # Validation
    validatetable = valdat[ids2,]
    
    output2 = matrix(rep(0,3),1,3)
    colnames(output2) = c('R2 Adjusted','Regression Coef','P-value')
    tem = data.frame(id=rownames(preds),prs = preds[,paste0('prs',Indx[trait,race])])
    tem$id = as.character(tem$id)
    rownames(tem) = tem$id
    tem = tem[validatetable$id,]
    
    prstable = merge(tem,validatetable,by='id')
    prstable = prstable[complete.cases(prstable),]
    if ((nrow(prstable)>0)&(sum(prstable$prs) != 0)){
      formula.prs=formula(paste0(paste('y', paste(c('age','gender',paste0('pc',1:10)),collapse="+"), sep='~')))
      fit = lm(formula.prs, data=prstable)
      prstable$res = fit$residuals
      prstable$pred.prs = predict(fit.train, data.frame(prs = prstable$prs))
      fit = lm(res~pred.prs, data=prstable)
      output2[1,'Regression Coef'] = coefficients(fit)['pred.prs']
      output2[1,'R2 Adjusted'] = summary(fit)$r.squared # (coefficients(fit)['prs'])^2/var(na.omit(prstable[,'res']))
      output2[1,'P-value'] = summary(fit)$coefficients['pred.prs','Pr(>|t|)']
    }
    out[1,paste(c('R2 Adjusted','Regression Coef','P-value'),race)] = output2[1,c('R2 Adjusted','Regression Coef','P-value')]
    
    Out[trait,paste('R2',race)] = mean(out[,paste(c('R2 Adjusted'),race)])
    Out[trait,paste('SD',race)] = sqrt(var(out[,paste(c('R2 Adjusted'),race)]))
    Out[trait,paste('P-value',race)] = mean(out[,paste(c('P-value'),race)])
    print(Out[trait,])
    print(trait)
    pars[paste(trait,race),] = sets[Indx[trait,race],]
    R2[trait,race] = output2[1,'R2 Adjusted']
    save(pars,file=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/r2/tuned-pars-ldpred2-",race,"-",trait,".RData"))
  }
}

####### -------------------------------  tuned parameter:
save(R2, Out, Indx, file=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision1/ldpred2_sequence/r2/r2-ldpred2.RData"))
