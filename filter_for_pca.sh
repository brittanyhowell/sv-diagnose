#!/bin/bash
## Usage: Filter VCF for basic quality cutoffs for PCA 
## Author: Brittany Howell (bh10@sanger.ac.uk)
## Date: 6th May 2019

wholeVCF="/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/pca/BATCHmerge.cn.vcf.gz"

cd /lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/pca/filter

## Here I'll use bcftools filter to remove variants based on their AC
echo "/nfs/team151/software/bcftools_2018/bcftools-1.9/bcftools filter -e 'INFO/AC < 2' -O z -o ${wholeVCF%cn.vcf.gz}AC2.vcf.gz ${wholeVCF} "  > run_exclude.sh 
bsub -G team151 -o "%J_bcf_exclude.o" -R"span[hosts=1]" -R"select[mem>500] rusage[mem=500]" -M500 "bash run_exclude.sh"