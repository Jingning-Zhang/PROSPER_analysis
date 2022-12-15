

## update geno snp id
#library(bigreadr)
#library(dplyr)
#library(readr)
#prs_cs_ref = read.table("/dcs04/nilanjan/data/jzhang2/TOOLS/prscsx/LDref/snpinfo_mult_1kg_hm3",header=T)


bim <- fread2("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/EAS/tv_allchr.bim")
for(chr in 1:22){
  tmp <- inner_join(bim[bim$V1==chr,c(2,4)],prs_cs_ref[prs_cs_ref$CHR==chr,2:3],by=c("V4"="BP"))[,c("V2","SNP")]
  if(chr == 1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
}
write_tsv(res, "/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/tmp/EAS_update_rsid", col_names=F)

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EAS")
system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 1 ",
" --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/EAS/tv_allchr ",
" --keep /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/EAS/validation.id.txt ",
" --update-name /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/tmp/EAS_update_rsid ",
" --make-bed ",
" --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EAS/tuning_geno"))

system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 1 ",
" --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/EAS/tv_allchr ",
" --keep /dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/EAS/test.id.txt ",
" --update-name /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/tmp/EAS_update_rsid ",
" --make-bed ",
" --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/EAS/testing_geno"))
