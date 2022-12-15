#!/usr/bin/env bash
#$ -N clean_score
#$ -cwd
#$ -l mem_free=15G,h_vmem=15G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu

skipnode='!(compute-053|compute-054|compute-113)'

module load R/3.6.1

R CMD BATCH --no-save --no-restore 1_clean_score.R

## -hold_jid tv_*