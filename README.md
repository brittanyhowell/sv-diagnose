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
for sample in `ls *.vcf | sed 's/.gt.vcf//g'`
    do
        cat ${sample}.gt.vcf | grep -v "#"  | grep  -v "]"  | grep  -v "\["  | cut -f1,2,5,6,10 | cut -f1 -d: | sed 's/<//g' | sed 's/>//g' > info.A
        cat ${sample}.gt.vcf | grep -v "#" | grep  -v "]"  | grep  -v "\["  | cut -f8 | cut -f2 -d\; | sed 's/SVLEN=//g'  > info.B
        paste info.A info.B > /lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/BATCH1/${sample}_vars.txt
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