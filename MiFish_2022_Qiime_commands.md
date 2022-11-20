MiFish 2022 samples
================
Gemma Clucas
2022-11-15

## 1. Import the data into Qiime2

First, load qiime environment and cd to correct directory in the
terminal.

    cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/Matinicus-Rock-2021-Atlantic-Puffins
    conda activate qiime2-2021.4

The raw reads are saved on my solid state hardrive, Data_SS1. There are
nearly three whole plates, although the third plate has a few Troppy
samples on it.

    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/MiFish/ATPU_WillK2022_Plate1/reads \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path MiFish_2022/demux_plate1.qza
      
    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/MiFish/ATPU_WillK2022_Plate2/reads \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path MiFish_2022/demux_plate2.qza
      
    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/MiFish/ATPU_WillK2022_Plate3/reads \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path MiFish_2022/demux_plate3.qza  

Summarise read quality and number of reads for each plate.

    cd MiFish_2022/

    for K in {1..3}; do
      qiime demux summarize \
        --i-data demux_Plate$K.qza \
        --o-visualization demux_Plate$K.qzv
    done

Plates 2 and 3 had more reads and higher quality than plate 1. Plate 1
still looks ok though I think.

## 2. Trim primers using cutadapt

The MiFish sequences are:

F primer: GTCGGTAAAACTCGTGCCAGC (21 bp)  
R primer: CATAGTGGGGTATCTAATCCCAGTTTG (27 bp)

### Trim 3’ ends first

At the 3’ end of the read, the primer will have been read through after
reading the MiFish region. I need to be looking for the reverse
complement of the reverse primer in read 1 (—p-adapter-f) and the
reverse complement of the forward primer in R2 (—p-adapter-r).

F primer reverse complement: GCTGGCACGAGTTTTACCGAC  
R primer reverse complement: CAAACTGGGATTAGATACCCCACTATG

    for K in {1..3}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux_Plate$K.qza \
        --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
        --p-adapter-r GCTGGCACGAGTTTTACCGAC \
        --o-trimmed-sequences trimd_Plate$K.qza \
        --verbose > cutadapt_out_Plate$3.txt
    done

I accidentally messed up the verbose bit of this command, so it only
printed the output for the last plate (the others would have been
over-written). But this isn’t really a big deal.

Make new visualisations to see how many sequences are left and their
quality scores.

    for K in {1..3}; do
      qiime demux summarize \
        --i-data trimd_Plate$K.qza \
        --o-visualization trimd_Plate$K.qzv
    done 

All looks good.

### Trim 5’ ends of reads

All R1 should begin with the forward primer: GTCGGTAAAACTCGTGCCAGC (21
bases). All R2 should begin with the reverse primer:
CATAGTGGGGTATCTAATCCCAGTTTG (27 bases).

Trim these with the following commands:

    for K in {1..3}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences trimd_Plate$K.qza \
        --p-front-f GTCGGTAAAACTCGTGCCAGC \
        --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
        --o-trimmed-sequences trimd2_Plate$K.qza \
        --verbose > cutadapt_out2_Plate$K.txt
    done

To see how much data passed the filter for each sample:

    grep "Total written (filtered):" cutadapt_out2_Plate1.txt 
    grep "Total written (filtered):" cutadapt_out2_Plate2.txt
    grep "Total written (filtered):" cutadapt_out2_Plate3.txt

About 90% passed the filter here for all samples. That’s pretty normal.

## 3. Denoise with dada2

I am going to use the same settings that I found worked best for the
2017-2019 tern fecal samples, except I am adding the –p-min-overlap
parameter. If I don’t use this, I seem to be getting a load of rubbish
reads which are 250bp long and start with long strings of Cs. This
option is only available in qiime2-2021.4 and later. I didn’t get these
rubbish reads before, so I’m not sure what has changed, but the overlap
filter seems to fix it.

Note, this step can be a little slow to run.

    for K in {1..3}; do
      qiime dada2 denoise-paired \
        --i-demultiplexed-seqs trimd2_Plate$K.qza \
        --p-trunc-len-f 133 \
        --p-trunc-len-r 138 \
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

## 4. Merge across plates

Note that this requires metadata for each sample, contained in
metadata.txt.

    qiime feature-table merge \
      --i-tables table_Plate1.qza \
      --i-tables table_Plate2.qza \
      --i-tables table_Plate3.qza \
      --p-overlap-method sum \
      --o-merged-table merged-table.qza
      
    qiime feature-table summarize \
        --i-table merged-table.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization merged-table
        
    qiime feature-table merge-seqs \
      --i-data rep-seqs_Plate1.qza \
      --i-data rep-seqs_Plate2.qza \
      --i-data rep-seqs_Plate3.qza \
      --o-merged-data merged_rep-seqs.qza
      
    qiime feature-table tabulate-seqs \
      --i-data merged_rep-seqs.qza \
      --o-visualization merged_rep-seqs.qzv

## 5. Assign taxonomy using same database as for 2021 samples

Note, I am using the exact same database files that I used for the 2021
samples.

    conda activate qiime2-2019.4

    ../MiFish/mktaxa.py \
      ../MiFish/editing_database/ncbi-refseqs-withHuman.qza \
      ../MiFish/editing_database/ncbi-taxonomy-withHuman.qza \
      merged_rep-seqs.qza

Make the visualisation.

    qiime metadata tabulate \
      --m-input-file superblast_taxonomy.qza \
      --o-visualization superblast_taxonomy

## 6. Make some barplots

I want to see how many puffin, human, and unassigned sequences there
were.

    qiime taxa barplot \
      --i-table merged-table.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_before_filtering.qzv

As with the 2021 samples, quite a lot of the samples had 100% avian DNA,
which is annoying since we can’t recover any fish DNA from them. I
haven’t seen this in terns before. Perhaps the primers are a better
match to puffins? Or perhaps the birds hadn’t fed on anything recently?
Not sure.

## 7. Remove non-food reads

Filter out any sequences from the bird, mammals (human), and unnassigned
sequences since we’re not interested in these.

    qiime taxa filter-table \
      --i-table merged-table.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --p-exclude Unassigned,Aves,Mammalia \
      --o-filtered-table merged_table_noBirdsMammalsUnassigned.qza
      
    qiime feature-table summarize \
        --i-table merged_table_noBirdsMammalsUnassigned.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization merged_table_noBirdsMammalsUnassigned

Barplot having removed bird/human/unassigned DNA:

    qiime taxa barplot \
      --i-table merged_table_noBirdsMammalsUnassigned.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_noBirdsMammalsUnassigned.qzv

## 9. Calculate alpha rarefaction curves

Last year we found a depth of 600 was sufficient for recovering all of
the diversity in the samples. I should use this again this year to be
consistent, but I want to check this too.

Also, I get a seg fault in the older version, so reload qiime2-2021.4
first.

    conda activate qiime2-2021.4

    qiime taxa collapse \
      --i-table merged_table_noBirdsMammalsUnassigned.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --p-level 7 \
      --o-collapsed-table merged_table_noBirdsMammalsUnassigned_collapsed.qza

    qiime diversity alpha-rarefaction \
      --i-table merged_table_noBirdsMammalsUnassigned_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 100 \
      --p-max-depth 20000 \
      --o-visualization alpha-rarefaction-100-20000

The curve for the observed OTUs (species) increases from 100 to 2000,
and again above 9000, which is very high - maybe we’re picking up super
low concentrations of DNA from past meals at this depth.

To compare to last year, repeat the rarefaction for values less than
2000.

    qiime diversity alpha-rarefaction \
      --i-table merged_table_noBirdsMammalsUnassigned_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 100 \
      --p-max-depth 2500 \
      --o-visualization alpha-rarefaction-100-2500

Shannon diversity is still flat for the samples. The observed OTUs has
plateued by 400. This suggests that 600 that we used last year is fine
for now.

## 10. Rarefy to a deph of 600 and redo barplots

Rarefying all samples to a depth of 600 to be consistent with 2021
samples.

Note, this is done on the un-collapsed table.

    qiime feature-table rarefy \
      --i-table merged_table_noBirdsMammalsUnassigned.qza \
      --p-sampling-depth 600 \
      --o-rarefied-table merged_table_noBirdsMammalsUnassigned_rarefied600
      
    qiime taxa barplot \
      --i-table merged_table_noBirdsMammalsUnassigned_rarefied600.qza \
      --i-taxonomy superblast_taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_oBirdsMammalsUnassigned_rarefied600.qzv

I downloaded the CSV and sent it to Will.

I tried briefly to compare the 25 samples we had done last year and
again on the robot, but realised that the naming convention is not the
same, so ATPU1 from last year is not the same sample as ATPU2021_1. Will
look at this when there’s more time to figure out exactly which sample
is which.
