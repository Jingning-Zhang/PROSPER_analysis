#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -t 1-5
#$ -m e
#$ -M jzhan218@jhu.edu


cp without_UKB_HDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AFR/HDL.gz
cp without_UKB_nonHDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AFR/nonHDL.gz
cp without_UKB_LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AFR/LDL.gz
cp without_UKB_TC_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AFR/TC.gz
cp without_UKB_logTG_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AFR/logTG.gz


cp HDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AMR/HDL.gz
cp nonHDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AMR/nonHDL.gz
cp LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AMR/LDL.gz
cp TC_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AMR/TC.gz
cp logTG_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/AMR/logTG.gz

cp without_UKB_HDL_INV_SAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/SAS/HDL.gz
cp without_UKB_nonHDL_INV_SAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/SAS/nonHDL.gz
cp without_UKB_LDL_INV_SAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/SAS/LDL.gz
cp without_UKB_TC_INV_SAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/SAS/TC.gz
cp without_UKB_logTG_INV_SAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/SAS/logTG.gz

cp without_UKB_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EUR/HDL.gz
cp without_UKB_nonHDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EUR/nonHDL.gz
cp without_UKB_LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EUR/LDL.gz
cp without_UKB_TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EUR/TC.gz
cp without_UKB_logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EUR/logTG.gz


cp HDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EAS/HDL.gz
cp nonHDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EAS/nonHDL.gz
cp LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EAS/LDL.gz
cp TC_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EAS/TC.gz
cp logTG_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz /dcs04/nilanjan/data/jzhang2/GLGC_cleaned/GLGC_raw/EAS/logTG.gz

