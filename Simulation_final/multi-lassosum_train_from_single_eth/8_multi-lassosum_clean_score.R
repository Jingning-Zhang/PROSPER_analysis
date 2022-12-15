# /dcs04/nilanjan/data/jzhang2/MEPRS/simcodes
##############################

rm(list=ls())

library(readr)
library(bigreadr)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

NCORES <- 3

#############################################################
#############################################################

M <- length(ethnic)
ethnic1 <- paste(ethnic,collapse = "_")

if(fixEUR=="T"){
  path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/4_multi-lassosum_by_chr_",para,"_fixEUR")
}else{
  path=paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum_train_from_single_eth/Setting_",setting,"/",ethnic1,"/rep_",rep,"/2_multi-lassosum_by_chr_",para)
}

for (mm in 1:M){

  for (chr in 1:22){

    snps_scale <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/snps_scale_rho_",rho,"_size_",size,"_chr",chr,".rds"))

    print(chr)

    # indx, indx_block, M, snp_list, Nsnps, ethnic, delta_best, summ_max, N,
    load(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/multi-lassosum/Setting_",setting,"/",ethnic1,"/rep_",rep,"/1_standard_data/otherinfo_rho_",rho,"_size_",size,"_chr",chr,".RData"))

    # ethnic, M, delta_best, runtime, res,
    load(paste0(path,"/tuning/rho_",rho,"_size_",size,"_chr",chr,".RData"))
    delta <- delta_best

    if(chr==1){
      alltuning <- length(res$c)
      #conv <- matrix(nrow=alltuning,ncol = 22)
      ## param
      param <- matrix(nrow=alltuning, ncol = (2*M+1+22))
      for (m in 1:M){ param[,m] <- delta[m] }
      param[,(M+1):(2*M)] <- t(res$lambda)
      param[,(2*M+1)] <- unlist(lapply(res$c, FUN = function (x){x[1,2]}))
    }

    ## snpsinfo
    snps <- unlist(snp_list)
    bim.ref <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",ethnic[mm],"/ref_chr", chr,".bim")) # !! may need to make sure the ref allele
    bim.ref <- bim.ref[match(snps, bim.ref$V2),]
    A1 <- bim.ref$V5
    m_snps <- !is.na(A1)

    snpsinfo <- data.frame(rsid=snps[m_snps], a1=A1[m_snps], stringsAsFactors = F)

    if( (fixEUR == "T") & (ethnic[mm]=="EUR") ){
      Beta_all <- matrix(nrow = sum(m_snps), ncol = 1)
      i <- 1
      Beta_all[,i] <- unlist(lapply(res$b[[i]], FUN=function (x) {x[,mm]}))[m_snps] * unlist(lapply(snps_scale, FUN=function (x) {x[,mm]}))[m_snps]
    }else{
      Beta_all <- matrix(nrow = sum(m_snps), ncol = alltuning)
      for (i in 1:alltuning){ Beta_all[,i] <- unlist(lapply(res$b[[i]], FUN=function (x) {x[,mm]}))[m_snps] * unlist(lapply(snps_scale, FUN=function (x) {x[,mm]}))[m_snps] } # dosage effect
    }

    param[,(2*M+1+chr)] <- apply(Beta_all, MARGIN = 2, FUN = function (x){mean(x!=0)})
    colnames(param) <- c(paste0("delta",1:M), paste0("lambda",1:M), "c", paste0("sparsity_chr",1:22))

    Beta_all[is.na(Beta_all)] <- 0
    Beta_all[Beta_all > 10] <- 0; Beta_all[Beta_all < -10] <- 0

    save(Beta_all, snpsinfo,
         file = paste0(path,"/Results/clean_results/beta+snps_rho_",rho,"_size_",size,"_",ethnic[mm],"_chr_",chr,".RData"))

    if(chr==22){
      write_tsv(as.data.frame(param), paste0(path,"/Results/clean_results/para_rho_",rho,"_size_",size,"_",ethnic[mm],".txt"))
    }

    rm(list = c("snps","m_snps","bim.ref","A1","snpsinfo","Beta_all","res","snps_scale"))
  }

  for (chr in 1:22){

    print(chr)
    load(paste0(path,"/Results/clean_results/beta+snps_rho_",rho,"_size_",size,"_",ethnic[mm],"_chr_",chr,".RData"))

    if((fixEUR == "T") & (ethnic[mm]=="EUR")){ Beta_all <- Beta_all[,1,drop=F] }
    Beta_all <- signif(Beta_all, digits = 3) # just save 3 significant

    tmp <- apply(Beta_all, MARGIN=1, function(x){sum(x!=0)}); m <- !(tmp==0)
    Beta_all <- Beta_all[m,,drop=F]
    snpsinfo <- snpsinfo[m,,drop=F]

    Beta_all <- data.frame(snpsinfo, Beta_all)
    fwrite2(Beta_all,
            paste0(path,"/Results/clean_results/Beta_rho_",rho,"_size_",size,"_",ethnic[mm],"_chr_",chr,".txt"),
            col.names = F, sep="\t", nThread=NCORES)
    rm(list = c("Beta_all","tmp","snpsinfo"))

  }

}
