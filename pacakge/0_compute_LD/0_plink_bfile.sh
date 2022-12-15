#!/usr/bin/env bash
#$ -N sumdata
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jzhan218@jhu.edu


mkdir /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno/EAS/

for (( i = 1; i < 23; i++ )); do
/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir/plink2 \
  --bfile /dcs04/nilanjan/data/jzhang2/1000G/GRCh37/EAS/chr${i} \
  --extract /dcl01/chatterj/data/jin/MEGA/mega-hm3-rsid.txt \
  --make-bed \
  --out /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/ref_geno/EAS/ref_chr${i} \
  --threads 1
done

/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir/plink2 \
  --bfile /dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/prepare_data/LD_block/AFR/tmp/byblock/chr1/chr1_247344518_249239466 \
  --freq \
  --out /users/jzhang2/tmp_tmp/freq