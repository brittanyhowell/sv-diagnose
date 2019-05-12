# Scripts written for the purpose of figuring out why my plots look so strange

## Where we are

I'm going to work in here: ```/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/```
and I'll use an interactive session: ```bsub -G team151 -Is -M 6000 -R'select[mem>6000] rusage[mem=6000]' -q yesterday  R-3.4.0```

## Step one: Plot the number of variants, but not as a summary, as a set of diagnostic plots

It's important I have as many  metrics as possible to plot. But I'm not entirely sure how to include all of that as a single data frame. Also I don't think it's entirely smart to do so.

### Plot one: Plot the number of variants per person, colour according to batch

Data required: variant info, genotypes, list of which samples were in which batches.

I already have genotypes per person, I need to count how many of each geno per person and save it somewhere logical.

#### How many variants/genotypes per person - PostMerge

First I'll bring over a file containing genotypes per sample

``` bash
ln -s /lustre/scratch119/humgen/projects/cnv_15x/svtools/Phase-I/BATCHmerge/WDL_scripts/cromwell-executions/Post_Merge_SV/48e5155e-4aeb-4718-b698-c86a256855ac/call-Sort_Index_VCF/execution/BATCHmerge.AF.CN.GT.nosq.txt ./BATCHmerge.AF.CN.GT.nosq.txt
```

Now I'll start with an adapted version of analyseVCF_per_variant.R - plot_number_per_sample.R

**Output file:** ```frequency_genotypes.txt```

Now I have to bind this to which batch the variants belong to.

#### In which batch do the variants belong

Ugh I have deleted the individual filelists per batch. Nevermind I'll make them again: 

``` bash
cd /lustre/scratch119/humgen/projects/cnv_15x/svtools/Phase-I/PostMerge_vcf/
head -n4000 BATCH1.vcf | grep "SAMPLE=" | sed 's/##SAMPLE=<ID=//g'| sed 's/>//g' | awk '{print $0 "\tBATCH1"}' > sampleList_batch1.list
head -n4000 BATCH2.vcf | grep "SAMPLE=" | sed 's/##SAMPLE=<ID=//g'| sed 's/>//g' | awk '{print $0 "\tBATCH2"}' > sampleList_batch2.list
head -n4000 BATCH3.vcf | grep "SAMPLE=" | sed 's/##SAMPLE=<ID=//g'| sed 's/>//g' | awk '{print $0 "\tBATCH3"}' > sampleList_batch3.list
head -n4000 BATCH4.vcf | grep "SAMPLE=" | sed 's/##SAMPLE=<ID=//g'| sed 's/>//g' > tmp.tmp
cat tmp.tmp missing | awk '{print $0 "\tBATCH4"}' > sampleList_batch4.list
rm tmp.tmp
head -n4000 BATCH5.vcf | grep "SAMPLE=" | sed 's/##SAMPLE=<ID=//g'| sed 's/>//g' | awk '{print $0 "\tBATCH5"}' > sampleList_batch5.list
cat sampleList_batch*list > sampleList_batchname.txt
mv sampleList_batchname.txt /lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/
```

Now I have a table to read into R --- Back to plot_number_per_sample.R

### How many variants/genotypes per person - PreMerge

I have worked out how many variants per sample after regenotyping and it is wacky. So I am going to return to the original VCFs, and count again. The following loop is run inside of a directory containing all of the VCFs.

#### read a vcf

``` bash
cd /lustre/scratch119/humgen/projects/cnv_15x/svtools/Phase-I/PreMerge_vcf/BATCH5
for sample in `ls *.vcf | sed 's/.gt.vcf//g'`
    do
        cat ${sample}.gt.vcf | grep -v "#"  | grep  -v "]"  | grep  -v "\["  | cut -f1,2,5,6,10 | cut -f1 -d: | sed 's/<//g' | sed 's/>//g' > info.A
        cat ${sample}.gt.vcf | grep -v "#" | grep  -v "]"  | grep  -v "\["  | cut -f8 | cut -f2 -d\; | sed 's/SVLEN=//g'  > info.B
        paste info.A info.B > /lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/BATCH5/${sample}_vars.txt
    done
```

The resulting text files look like this:

``` bash
chr1    1657106    DUP    52.68    0/1    65269
chr1    2652330    DEL    0.00     0/0    -2178
chr1    2122016    DEL    251.80   0/1    -2208
chr1    3203223    DEL    2.88     0/0    -360
```

Which I will now feed back into plot_number_per_sample.R

### Results from PreMerge

Well, they certainly looked more normal.
Now what.
Something that might be having an effect is the number of variants in the merge file.

#### Number of non BND lines in the merge vcf

1. 1304385
2. 1275944
3. 2600882
4. 2851188
5. 2614693

Note this is the sorted file, not the merged file. I had to delete the merge file. Well, I didn't have to, but I did because I had to delete a lot of intermediate files in a rush to get back some free space on the disk.

#### regenerate them

If I could separate out the l-merge function from the MergeWDL, then I could specifically rerun the relevant commands.
I suppose this isn't so hard (see [rerun_lmerge.sh](https://github.com/brittanyhowell/sv-diagnose/))
I ran l-merge from the Merge.WDL on the sorted vcfs, and now have one merge file per batch.

#### Analyse

Reading:

``` bash
# Unzip all:
gzip -d ${batch}.gz

# Pull out some feats:
for batch in `ls *5*.vcf`
    do
        cat ${batch} | grep -v "#"  | grep  -v "]"  | grep  -v "\["  | cut -f1,2,5,6  | sed 's/<//g' | sed 's/>//g' > info.A
        cat ${batch} | grep -v "#" | grep  -v "]"  | grep  -v "\["  | cut -f8 | perl -pe 's/SVLEN=/\nSVLEN=/g' |  perl -pe 's/;CIPOS/\n/g' | perl -pe 's/;END/\n/g' | grep "SVLEN" |  perl -pe  's/SVLEN=//g'  > info.B
        paste info.A info.B > ${batch%.vcf}_vars.txt
    done
```

See Merge vcf characteristics. in plot_number_per_sample.R

## Checking the Pre_Post_Merge VCFs

Checking these will allow me to know how many per variant and also I can plot more metrics maybe

### What to extract

The VCF options:
```CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  EGAN00001345375```

With the Pre_Merge VCFs, I just read over it with this:

``` bash
# Make the output folder
mkdir -pv /lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/Pre_Post_Merge/BATCH3/

cd /lustre/scratch119/humgen/projects/cnv_15x/svtools/Phase-I/PrePostMerge_vcf/BATCH3
for sample in `ls *.vcf | sed 's/.gt.vcf//g'`
    do
        cat ${sample}.gt.vcf | grep -v "#"  | grep  -v "]"  | grep  -v "\["  | cut -f1,2,5,6,10 | cut -f1 -d: | sed 's/<//g' | sed 's/>//g' > info.A
        cat ${sample}.gt.vcf | grep -v "#" | grep  -v "]"  | grep  -v "\["  | cut -f8 | cut -f2 -d\; | sed 's/SVLEN=//g'  > info.B
        paste info.A info.B > /lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/Pre_Post_Merge/BATCH3/${sample}_vars.txt
    done

    bsub -o "makeSums_batch3_%J.o" "bash makeSumsB3.sh"
```

Which is fine, and what I'm going to start with, however, I could pull out any of THESE:

``` bash
bcftools query GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB
```

Ok really though once PPMs are all unzipped, I can read them into the old faithful R script

## PCA

### Preparing the data

I think the PCA should be split into two: stats per variant and stats per individual.
Ideally, I would have a large text file which contains all variants, from all

1. bcftools query each and every table to pull out relevant info fields
2. read into R for processing

Going to paste on a batch number and randomly select a few vars from each batch:

``` bash
cd /lustre/scratch119/humgen/projects/cnv_15x/svtools/Phase-I/PostMerge_vcf
bsub -o "zip_me_%J.o" "bash zipMe.sh"

bcftools view -Oz -o  BATCH1.vcf.gz BATCH1.vcf
bcftools view -Oz -o  BATCH2.vcf.gz BATCH2.vcf
bcftools view -Oz -o  BATCH3.vcf.gz BATCH3.vcf
bcftools view -Oz -o  BATCH4.vcf.gz BATCH4.vcf
bcftools view -Oz -o  BATCH5.vcf.gz BATCH5.vcf

bcftools query  -f '%CHROM\t%POS\t%ID\t%QUAL\t%SVTYPE\t%AF\t%NSAMP\t%MSQ\n' BATCH5.vcf.gz > BATCH5_from_query.txt
cat BATCH5_from_query.txt | grep -v "BND" |awk '{print $0 "\t BATCH5"}' > BATCH5.txt

shuf -n 1000 BATCH5.txt > BATCH5_rand.txt

cat BATCH1_rand.txt BATCH2_rand.txt BATCH3_rand.txt BATCH4_rand.txt BATCH5_rand.txt > BATCHall_rand.txt

```

## A different type of PCA

Since I can't really figure out a nice way to come up with my own eigenvectors (fairly sure this is wrong _anyway_) I will just give everything to plink and throw everything to the wind.

### Make Plink format, make PCA

I used a variation of this code, with differences in ```maf```, ```hwe``` and ```geno``` for different results.

``` bash
bsub -G team151 -Is -M 6000 -R'select[mem>6000] rusage[mem=6000]' -q yesterday bash


/lustre/scratch118/infgen/team133/db22/software/plink2/plink2 --vcf BATCHmerge.GT.vcf.gz --out reduce --make-bed

# /lustre/scratch118/infgen/team133/db22/software/plink2/plink2 --vcf BATCHmerge.GT.vcf.gz  --make-bed --out LD/LD

# /lustre/scratch118/infgen/team133/db22/software/plink2/plink2 --bfile LD/LD --indep-pairwise 50 10 0.15 --out LD/LD

# /lustre/scratch118/infgen/team133/db22/software/plink2/plink2 --bfile LD/LD --extract LD/LD.prune.in --out LD/LDout --make-bed


/lustre/scratch118/infgen/team133/db22/software/plink2/plink2 --bfile reduce --pca --out reduce_pca
```

R code available in pca.R.
I did some basic filtering using Plink, as above, but I expanded on this in the script filter_for_pca.sh

### Removing many fields from the VCF

To see what is affecting the PCA, I'm going to strip back the VCF so it plots with fewer metrics

``` bash
/nfs/team151/software/bcftools_2018/bcftools-1.9/bcftools annotate -O z -o BATCHmerge.reduce.vcf.gz -x ^INFO/AF,^INFO/SVLEN,^INFO/MSQ,^FORMAT/GT BATCHmerge.cn.vcf.gz
```

### Metadata

I need to colour my samples by more than just batch. I'm going to work on a metadata table for them. - I stalled here because I didn't like the idea of summary statistics and counting them, and was running out of ideas.

#### Information needed

1. number of sites per sample - of each type
2. Average quality of samples
3. Number of bases called in that person.

## A kind of association

To test for systematic biases, I am going to conduct a case control analysis between batches, and use some diagnostic plots to see what to blame.

## Further action

1. Plot QUAL dist <2000
1. Merge batches 1&2, regenotype.
2. Make PCA of AF <%5, <1%, of all samples
2. Check the impact of the SU filter

1. Look at new code changes
2. Test new code with downloaded batch
3. set up regenotyping for batches 1&2

4. Make a metadata table