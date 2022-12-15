#!/usr/bin/env bash
#$ -N GRCh37
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

module load conda_R/4.0

R CMD BATCH --no-save --no-restore 0_clean_mage_MAF_GRCh37.R
