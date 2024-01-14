
library(bigreadr)
library(readr)
library(bigsnpr)

parameters = expand.grid(race = c('AFR','AMR','EAS','EUR','SAS'),
                         setting = 1:5, rho = 1:3, size = 1:4, rep = 1)

args <- commandArgs(T)
for(i in 1:length(args)){ eval(parse(text=args[[i]])) }

race = as.character(parameters[as.numeric(temp1),1])
setting = as.integer(parameters[as.numeric(temp1),2])
rho = as.integer(parameters[as.numeric(temp1),3])
size = as.integer(parameters[as.numeric(temp1),4])
rep = as.integer(parameters[as.numeric(temp1),5])


print(setting)
print(rho)
print(size)
print(rep)

print(race)

NCORES <- 1

sumstats <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/summary_stat_bigsnpr/",race,"/rho_",rho,"_size_",size,"_rep_",rep,"_GA_",setting,".txt"))

# -------- QC of summary statistics and reference data --------

# Match SNPs between summary statistics and reference data
ref <- fread2(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/",race,"/ref_allchr.bim"))

snp <- intersect(ref$V2,sumstats$rsid)
ref <- ref[match(snp,ref$V2),]
sumstats <- sumstats[match(snp,sumstats$rsid),]
sumstats$pos <- ref$V4

flip <- ref$V5 != sumstats$a1
if(sum(flip)!=0){
  sumstats$beta[flip] <- -sumstats$beta[flip]
  sumstats$a1[flip] <- ref$V5[flip]
  sumstats$a0[flip] <- ref$V6[flip]
}
rm(list = c("ref"))

map <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation_final/LD_bigsnpr/",race,"/map.rds"))
info_snp <- snp_match(sumstats, map)

rm(list=c("sumstats"))

info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp))
#sd_ldref <- with(info_snp, sqrt(2 * a1_af * (1 - a1_af)))
##sd_ss <- with(info_snp, sqrt(2 * a1_sumdata_af * (1 - a1_sumdata_af)))
#sd_ss <- with(info_snp, 2/sqrt(n_eff*beta_se^2+beta^2))
##is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
#is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss < 0.1 | sd_ldref < 0.1
#df_beta <- info_snp[!is_bad, ]
df_beta <- info_snp

dir.create( paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",race,"/rep_",rep,"/summary_stat"), recursive=T )
write_tsv(df_beta, paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/revision2/update_bigsnpr_results/Simulation/Setting_",setting,"/",race,"/rep_",rep,"/summary_stat/rho_",rho,"_size_",size,".txt"))

