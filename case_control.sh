#!/bin/bash
## Usage: Compare batches as if one was a case and one was a control
## Author: Brittany Howell (bh10@sanger.ac.uk)
## Date: 9th May 2019


VCFTOOLS=/nfs/team151/software/vcftools/bin/vcftools
DEBUGDIR=/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/
wkdir=${DEBUGDIR}/caseControl


BATCHA="5"
# BATCHB="5"
# keepFILE=${wkdir}/samples_keep_${BATCHA}_${BATCHB}.tmp
keepFILE=${wkdir}/samples_keep_${BATCHA}.tmp
# ln -s /lustre/scratch119/realdata/mdt2/projects/cnv_15x/svtools/Phase-I/BATCHmerge/WDL_scripts/cromwell-executions/Post_Merge_SV/48e5155e-4aeb-4718-b698-c86a256855ac/call-Sort_Index_VCF/execution/BATCHmerge.cn.vcf.gz  BATCHmerge.cn.vcf.gz 
wholeVCF="BATCHmerge.cn.vcf.gz"
outFILEPref="BATCHmerge.${BATCHA}" 
# outFILEPref="BATCHmerge.${BATCHA}.${BATCHB}" 



cd ${wkdir}
    
    # First, assemble VCF with only samples from useful batches.
        fileListBATCHA=${DEBUGDIR}/sampleList_batch${BATCHA}.list
        # fileListBATCHB=${DEBUGDIR}/sampleList_batch${BATCHB}.list

        cat ${fileListBATCHA} ${fileListBATCHB} | cut -f1 > ${keepFILE}
        cat ${fileListBATCHA}  | cut -f1 > ${keepFILE}

    # Run vcftools to get rid of unwanted samples
       echo "/nfs/team151/software/bcftools_2018/bcftools-1.9/bcftools view -S ${keepFILE} -O z -o ${wkdir}/${outFILEPref}   ${wholeVCF}" > filter_samples_${BATCHA}.sh #filter_samples_${BATCHA}_${BATCHB}.sh

        bsub -o filter_samples_%J.o "bash filter_samples_${BATCHA}.sh"

    # Convert to Plink format
        /lustre/scratch118/infgen/team133/db22/software/plink2/plink2 --vcf BATCHmerge.1.2.vcf.gz --out BATCHmerge.1.2 --make-bed

    # Change .fam file so A is cases and B is controls 
        mv BATCHmerge.${BATCHA}.${BATCHB}.fam BATCHmerge.${BATCHA}.${BATCHB}.real.fam
        touch BATCHmerge.${BATCHA}.${BATCHB}.fam
        
        # Controls (1)- BATCHA
        for sample in `cat ${fileListBATCHA} | cut -f1 | head -n200`
        do 
            awk -v sample="$sample" '{ if ($2 == sample) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t1" } }' BATCHmerge.${BATCHA}.${BATCHB}.backup >> BATCHmerge.${BATCHA}.${BATCHB}.fam
        done
        
        # Cases (2) - BATCHB
        for sample in `cat ${fileListBATCHA} | cut -f1 | tail -n200`
        do 
            awk -v sample="$sample" '{ if ($2 == sample) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t2" } }' BATCHmerge.${BATCHA}.${BATCHB}.backup >> BATCHmerge.${BATCHA}.${BATCHB}.fam
        done
    
    # run association
        /lustre/scratch118/infgen/team133/db22/software/plink1.9/plink --bfile BATCHmerge.1.2 --out BATCHmerge.1 --allow-no-sex  --assoc --adjust qq-plot




    
    