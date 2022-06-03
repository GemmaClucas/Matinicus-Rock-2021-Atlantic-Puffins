18S ATPU analysis in Qiime2
================
Gemma Clucas
6/3/2022

## 1\. Import the data into Qiime2

First, load qiime environment and cd to correct directory in the
terminal.

    cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/Matinicus-Rock-2021-Atlantic-Puffins
    conda activate qiime2-2021.4

The raw reads are saved on my solid state hardrive, Data\_SS1. There are
two whole plates, plates 1 and 2, and then some samples on
Fecal\_Plate\_15. I moved these extra files into their own directory so
as not to confuse them with the terns on that plate.

``` 
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/Data_SS1/18S/ATPU_WillK_2021/2021_ATPU_Plate1_WillK/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path 18S/demux_plate1.qza
  
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/Data_SS1/18S/ATPU_WillK_2021/2021_ATPU_Plate2_WillK/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path 18S/demux_plate2.qza
  
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/Data_SS1/18S/ATPU_WillK_2021/2021_ATPU_WillK_Plate15extras/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path 18S/demux_plate3.qza  
```

Move to the diectory where all the next steps will take place. Summarise
read quality and number of reads for each plate.

    cd 18S/
    
    for K in {1..3}; do
      qiime demux summarize \
        --i-data demux_plate$K.qza \
        --o-visualization demux_Plate$K.qzv
    done

View `.qzv` files at [view.qiime2.org](view.qiime2.org). The quality of
the reads looks good. Most samples have tens of thousands of reads,
which should be enough.

## 2\. Trim primers using cutadapt

The MiFish sequences are:

F primer: GGTCTGTGATGCCCTTAGATG (21 bp) R primer: GGTGTGTACAAAGGGCAGGG
(20 bp)

### Trim 3’ ends first

At the 3’ end of the read, the primer will have been read through after
reading the 18S region. I need to be looking for the reverse complement
of the reverse primer in read 1 (—p-adapter-f) and the reverse
complement of the forward primer in R2 (—p-adapter-r).

F primer reverse complement: CATCTAAGGGCATCACAGACC R primer reverse
complement: CCCTGCCCTTTGTACACACC

    for K in {1..3}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux_plate$K.qza \
        --p-adapter-f CCCTGCCCTTTGTACACACC \
        --p-adapter-r CATCTAAGGGCATCACAGACC \
        --o-trimmed-sequences trimd_Plate$K.qza \
        --verbose > cutadapt_out_Plate$K.txt
    done

To see how much data passed the filter for each sample:

    grep "Total written (filtered):" cutadapt_out_Plate1.txt 
    grep "Total written (filtered):" cutadapt_out_Plate2.txt
    grep "Total written (filtered):" cutadapt_out_Plate3.txt

About 76% of the data is passing the filters consistently across all
samples, which is a good sign.

Make new visualisations to see how many sequences are left and their
quality scores.

``` 
for K in {1..3}; do
  qiime demux summarize \
    --i-data trimd_Plate$K.qza \
    --o-visualization trimd_Plate$K.qzv
done 
```

Looking at the trimd.qzv files, things look good.

### Trim 5’ ends of reads

All R1 should begin with the forward primer: GGTCTGTGATGCCCTTAGATG (21
bases). All R2 should begin with the reverse primer:
GGTGTGTACAAAGGGCAGGG (20 bases).

Trim these with the following commands:

    for K in {1..3}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences trimd_Plate$K.qza \
        --p-front-f GGTCTGTGATGCCCTTAGATG \
        --p-front-r GGTGTGTACAAAGGGCAGGG \
        --o-trimmed-sequences trimd2_Plate$K.qza \
        --verbose > cutadapt_out2_Plate$K.txt
    done

Not going to make the qzv files for now.

To see how much data passed the filter for each sample:

    grep "Total written (filtered):" cutadapt_out2_Plate1.txt 
    grep "Total written (filtered):" cutadapt_out2_Plate2.txt
    grep "Total written (filtered):" cutadapt_out2_Plate3.txt

About 90% passed the filter here for all samples. That’s pretty normal.

## 3\. Denoise with dada2

With previous 18S data, truncating the length to 150 seemed to work
well. I haven’t used the –p-min-overlap setting before, as this is a new
addition to qiime2, so I am going to try with 50 (which is what I use
for the MiFish amplicon) and see how it goes.

Note, this step can be a little slow to run.

    for K in {1..3}; do
      qiime dada2 denoise-paired \
        --i-demultiplexed-seqs trimd2_Plate$K.qza \
        --p-trunc-len-f 150 \
        --p-trunc-len-r 150 \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-min-overlap 50 \
        --p-n-threads 16 \
        --o-representative-sequences rep-seqs_Plate$K \
        --o-table table_Plate$K \
        --o-denoising-stats denoise_Plate$K
    done

Create visualizations for the denoising stats.

    for K in {1..3}; do  
      qiime metadata tabulate\
        --m-input-file denoise_Plate$K.qza\
        --o-visualization denoise_Plate$K.qzv
    done

Seems like \>90% of the data is getting through all the denoising and
filtering done in this step, so that’s great.

To view the rep-seqs (only running this for one plate, will view all
after merging plates)

    qiime feature-table tabulate-seqs \
      --i-data rep-seqs_Plate1.qza \
      --o-visualization rep-seqs_Plate1

This looks good. All sequences are ~169-171bp. Using 50 as a minimum
overlap seems like a good setting to keep.

## 4\. Merge across plates

I need to merge both the feature tables, which contain the counts of
each feature, and the rep-seqs, which contain the actual sequence for
each feature.

Note that this requires metadata for each sample, contained in
metadata.txt. At the moment this is prety sparse as I didn’t have time
to add it all, but if you populate it with breeding stage, breeding
period, substrate etc then it’s easy to visalise those different groups
in the later steps of this analysis, so I would recommend editing the
metadata file. You can add any number of new columns to add info about
the samples. It just has to be saved as a tab-separated file.

    qiime feature-table merge \
      --i-tables table_Plate1.qza \
      --i-tables table_Plate2.qza \
      --i-tables table_extras.qza \
      --i-tables ../../2019_GoM_FecalMetabarcoding/Mifish/Puffin_tests/table_Puffin_tests.qza \
      --p-overlap-method sum \
      --o-merged-table merged-table.qza
      
    qiime feature-table summarize \
        --i-table merged-table.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization merged-table
        
    qiime feature-table merge-seqs \
      --i-data rep-seqs_Plate1.qza \
      --i-data rep-seqs_Plate2.qza \
      --i-data rep-seqs_extras.qza \
      --i-data ../../2019_GoM_FecalMetabarcoding/Mifish/Puffin_tests/rep-seqs_Puffin_tests.qza \
      --o-merged-data merged_rep-seqs.qza
      
    qiime feature-table tabulate-seqs \
      --i-data merged_rep-seqs.qza \
      --o-visualization merged_rep-seqs.qzv