

## update geno snp id
#library(bigreadr)
#library(dplyr)
#library(readr)
#prs_cs_ref = read.table("/dcs04/nilanjan/data/jzhang2/TOOLS/prscsx/LDref/snpinfo_mult_1kg_hm3",header=T)

a <- readLines("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/AFR/validation.id.txt")
a <- a[1:1000]
writeLines(a, "/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/example_in_github/id/AFR_tuning.id.txt")

a <- readLines("/dcs04/nilanjan/data/jzhang2/MEPRS/Simulation/geno/mega_matched_summary_stat/AFR/test.id.txt")
a <- a[1:1000]
writeLines(a, "/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/example_in_github/id/AFR_testing.id.txt")

dir.create("/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/example_in_github/sample_data/AFR")
system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 1 ",
" --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/tuning_geno ",
" --keep /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/example_in_github/id/AFR_tuning.id.txt ",
" --make-bed ",
" --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/example_in_github/sample_data/AFR/tuning_geno"))

system(paste0("/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2 --threads 1 ",
" --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/sample_data/AFR/testing_geno ",
" --keep /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/example_in_github/id/AFR_testing.id.txt ",
" --make-bed ",
" --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/example/example_in_github/sample_data/AFR/testing_geno"))
