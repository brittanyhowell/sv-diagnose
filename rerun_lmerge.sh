#!/bin/bash
## Usage: These commands were given to regenerate lmergeVCFs. 
##        They are here for reference, not intention to _reuse_
## Author: Brittany Howell (bh10@sanger.ac.uk)
## Date: 29th April 2019

    zcat BATCH1.lsort.vcf.gz\
      | svtools lmerge \
      -i /dev/stdin \
      -f 20 \
      | bgzip -c \
      > BATCH1.lmerge.vcf.gz

          zcat BATCH2.lsort.vcf.gz\
      | svtools lmerge \
      -i /dev/stdin \
      -f 20 \
      | bgzip -c \
      > BATCH2.lmerge.vcf.gz

          zcat BATCH3.lsort.vcf.gz\
      | svtools lmerge \
      -i /dev/stdin \
      -f 20 \
      | bgzip -c \
      > BATCH3.lmerge.vcf.gz

          zcat BATCH4.lsort.vcf.gz\
      | svtools lmerge \
      -i /dev/stdin \
      -f 20 \
      | bgzip -c \
      > BATCH4.lmerge.vcf.gz

          zcat BATCH5.lsort.vcf.gz\
      | svtools lmerge \
      -i /dev/stdin \
      -f 20 \
      | bgzip -c \
      > BATCH5.lmerge.vcf.gz

bsub -R"span[hosts=1]" -R"select[mem>10000] rusage[mem=10000]" -M10000 -o "mergeBatch1_%J.o" "bash b1merge.sh"
bsub -R"span[hosts=1]" -R"select[mem>10000] rusage[mem=10000]" -M10000 -o "mergeBatch2_%J.o" "bash b2merge.sh"
bsub -R"span[hosts=1]" -R"select[mem>10000] rusage[mem=10000]" -M10000 -o "mergeBatch3_%J.o" "bash b3merge.sh"
bsub -R"span[hosts=1]" -R"select[mem>10000] rusage[mem=10000]" -M10000 -o "mergeBatch4_%J.o" "bash b4merge.sh"
bsub -R"span[hosts=1]" -R"select[mem>10000] rusage[mem=10000]" -M10000 -o "mergeBatch5_%J.o" "bash b5merge.sh"




c("1304385", "1275944", "2600882", "2851188", "2614693")