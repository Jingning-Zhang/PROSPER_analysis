args = commandArgs(trailingOnly = T)
#i represent ethnic group
#j represent chromosome
#l represent the causal SNPs proportion
#m represent the training sample size
#i_rep represent simulation replication
#i1 represent the genetic architecture

i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
j = as.numeric(args[[3]])
v = as.numeric(args[[4]])
#j = as.numeric(args[[3]])
#m = as.numeric(args[[3]])
#i_rep = as.numeric(args[[4]])
#i1 = as.numeric(args[[5]])
library(data.table)
#install.packages("dplyr")
#install.packages("vctrs")
library(dplyr)

eth <- c("EUR","AFR","AMR","EAS","SAS")
trait <- c("any_cvd","depression",
           "heart_metabolic_disease_burden",
           "height",
           "iqb.sing_back_musical_note",
           "migraine_diagnosis",
           "morning_person")

sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
dir.create(paste0('/lscratch/',sid,'/test/1KGLD'),showWarnings = FALSE)
#copy the LD reference data to lscratch
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_eur ",temp.dir,"1KGLD"))
system(paste0("cp -r /data/zhangh24/software/PRScsx/1KGLD/ldblk_1kg_",tolower(eth[i])," ",temp.dir,"1KGLD"))
system(paste0("cp /data/zhangh24/software/PRScsx/1KGLD/snpinfo_mult_1kg_hm3 ",temp.dir,"1KGLD"))
sample_size <- as.data.frame(fread("/data/zhangh24/multi_ethnic/data/23_sample_size.csv",header=T))
setwd("/data/zhangh24/multi_ethnic/")
#temp.dir = paste0("/fastscratch/myscratch/hzhang1/ARIC/",trait[l],"/",eth[i],"/")
data.dir = "/data/zhangh24/multi_ethnic/data/cleaned/"
#out.dir = paste0("/data/zhangh24/multi_ethnic/result/cleaned/clumping_result/",method,"/",eth[i],"/",trait[l],"/")
#load EUR gwas summary statistics
sum.eur = as.data.frame(fread(paste0(data.dir,eth[1],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
# write.table(sum.data.MAF,file = paste0("/lscratch/",sid,"/test/",eth[i],"_summary_out_MAF_rho_",l,"_size_",m,"_rep_",i_rep,"_GA_",i1,".out")
#             ,col.names = T,row.names = F,quote=F) 
#prepare association file for plink
sum.eur.select = sum.eur %>% 
  rename(SNP=rsid) %>% 
  select(SNP,A1,A2,BETA,P) 
write.table(sum.eur.select,file = paste0(temp.dir,"EUR_sumstats.txt"),row.names = F,col.names = T,quote=F)
#load target ethnic group data
sum.tar = as.data.frame(fread(paste0(data.dir,eth[i],"/sumdat/",trait[l],"_passQC_noNA_matchinfo_matchMAFnoNA_common_mega+hapmap3_cleaned.txt"),header=T))
sum.tar.select = sum.tar %>% 
  rename(SNP=rsid) %>% 
  select(SNP,A1,A2,BETA,P) 
write.table(sum.tar.select,file = paste0(temp.dir,eth[i],"_sumstats.txt"),row.names = F,col.names = T,quote=F)
#create bim for prscsx
bim = sum.tar %>% 
  mutate(V3=0) %>% 
  select(CHR,rsid,V3,BP,A1,A2)
write.table(bim,file = paste0(temp.dir,eth[i],"_genotype.bim"),row.names = F,col.names = F,quote=F)
phi = c(1E+00,1E-02,1E-04,1E-6)


bintrait <- c("any_cvd","depression",
             
             "iqb.sing_back_musical_note",
             "migraine_diagnosis",
             "morning_person")
contriat = c("heart_metabolic_disease_burden",
             "height")

if(trait[l]%in%bintrait){
  n_eur = as.integer(4*sample_size %>% 
    filter(eth=="EUR"&
             Disease==trait[l]) %>% select(N_effect))
  n_tar = as.integer(4*sample_size %>% 
                       filter(eth==eth[i]&
                                Disease==trait[l]) %>% select(N_effect))
}else{
  n_eur = as.integer(sample_size %>% 
                       filter(eth=="EUR"&
                                Disease==trait[l]) %>% select(N_control))
  n_tar = as.integer(sample_size %>% 
                       filter(eth==eth[i]&
                                Disease==trait[l]) %>% select(N_effect))
}


path_to_ref = paste0(temp.dir,"1KGLD")
path_to_bim = paste0(temp.dir,eth[i],"_genotype")
path_to_sum = paste0(temp.dir)

out.dir.prs = paste0("/data/zhangh24/multi_ethnic/result/cleaned/prs/PRSCSx/",eth[i],"/",trait[l],"/")




system(paste0("python /data/zhangh24/software/PRScsx/PRScsx.py", 
              " --ref_dir=",path_to_ref,
              " --bim_prefix=",path_to_bim,
              " --sst_file=",path_to_sum,"EUR_sumstats.txt,",path_to_sum,eth[i],"_sumstats.txt",
              " --n_gwas=",n_eur,",",n_tar,
              " --pop=EUR,",eth[i],
              " --chrom=",j,
              " --phi=",phi[v],  
              " --out_dir=",out.dir.prs,
              " --out_name=sum"))




# out.file = paste0("/data/zhangh24/multi_ethnic/result/cleaned/organize_prs/PRSCSx/")
# for(i in 1:5){
#   system(paste0("mkdir ",out.file,eth[i],"/"))
#   for(l in 1:7){
#     system(paste0("mkdir ",out.file,eth[i],"/",trait[l]))  
#   }
# }
