
rm(list=ls())

library(bigreadr)
library(readr)

ancestry <- fread2("/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/ethnic_background_from_ukbb.txt")

## SAS (Indian, Pakistani, Bangladeshi)
SAS <- ancestry[ancestry$`f.21000.0.0` %in% c(3001,3002,3003),]
tmp <- format(SAS$f.eid, trim = T, scientific = F)
write_tsv(data.frame(FID=tmp, IID=tmp), "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_SAS.id")

## EAS (Chinese)
EAS <- ancestry[ancestry$`f.21000.0.0` %in% c(5),]
tmp <- format(EAS$f.eid, trim = T, scientific = F)
write_tsv(data.frame(FID=tmp, IID=tmp), "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_EAS.id")

## AFR (Caribbean, African)
AFR <- ancestry[ancestry$`f.21000.0.0` %in% c(4001,4002),]
tmp <- format(AFR$f.eid, trim = T, scientific = F)
write_tsv(data.frame(FID=tmp, IID=tmp), "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_AFR.id")

## smaples to exclude /dcl01/arking/data/UK_Biobank/static/Phenotype/Samples.to.Exclude/w17731_20220222.csv

## White (British, Irish, Any other white background) (can be removed from genetic ancestry prediction)
White <- ancestry[ancestry$`f.21000.0.0` %in% c(1001,1002,1003),]
tmp <- format(White$f.eid, trim = T, scientific = F)
write_tsv(data.frame(FID=tmp, IID=tmp), "/dcs04/nilanjan/data/jzhang2/UKBB/phenotype/ethnic_background/self_reported_White.id")

