# Scripts written for the purpose of figuring out why my plots look so strange

## Where we are

I'm going to work in here: ```/lustre/scratch119/humgen/projects/cnv_15x/svtools/debug/```
and I'll use an interactive session: ```bsub -G team151 -Is -M 6000 -R'select[mem>6000] rusage[mem=6000]' -q yesterday  R-3.4.0```

## Step one: Plot the number of variants, but not as a summary, as a set of diagnostic plots

It's important I have as many  metrics as possible to plot. But I'm not entirely sure how to include all of that as a single data frame. Also I don't think it's entirely smart to do so.

### Plot one: Plot the number of variants per person, colour according to batch

Data required: variant info, genotypes, list of which samples were in which batches.

I already have genotypes per person, I need to count how many of each geno per person and save it somewhere logical.

#### How many variants/genotypes per person

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
fileList=/lustre/scratch119/humgen/projects/cnv_15x/svtools/Phase-I/BATCHmerge/input_lists/sampleList.list  
files=( $(cat $fileList))
echo ${files[@]:0:465}     | tr " " "\n" | awk '{print $0 "\tBATCH1"}' > sampleList_batch1.list
echo ${files[@]:465:465}   | tr " " "\n" | awk '{print $0 "\tBATCH2"}' > sampleList_batch2.list
echo ${files[@]:930:930}   | tr " " "\n" | awk '{print $0 "\tBATCH3"}' > sampleList_batch3.list
echo ${files[@]:1860:932}  | tr " " "\n" | awk '{print $0 "\tBATCH4"}' > sampleList_batch4.list
echo ${files[@]:2792:932}  | tr " " "\n" | awk '{print $0 "\tBATCH5"}' > sampleList_batch5.list

cat *_batch* > sampleList_batchname.txt
```

Now I have a table to read into R --- Back to plot_number_per_sample.R


