# Matinicus-Rock-2021-Atlantic-Puffins

Started 15th Dec 2021

Fecal metabarcoding analysis for samples collected on Matinicus Rock by Will Kennerley, 2021.

Load qiime environment and cd to correct directory.
```
cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/Matinicus-Rock-2021-Atlantic-Puffins
conda activate qiime2-2021.4
```

## 1. Import the data into Qiime2

It's saved on my solid state hardrive, Data_SS1. There are two whole plates, plates 1 and 2, and then some samples on Fecal_Plate_15. 

There were hidden files in the folders, which I found using ```ls -alh```. Then I removed them using ```rm ._ATP*```. If you don't remove them then the commands won't run. 

```
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/Data_SS1/MiFish/ATPU_WillK2021_Plate1/211203_A01346_0039_BHTJ7TDRXY_16Mer120321_GC_ATPU1/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_plate1.qza
  
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/Data_SS1/MiFish/ATPU_WillK2021_Plate2/211203_A01346_0039_BHTJ7TDRXY_16Mer120321_GC_ATPU2/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_plate2.qza
  
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /Volumes/Data_SS1/MiFish/ATPU_WillK2021_extras/reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt\
  --output-path MiFish/demux_extras.qza  

```  
Move to the diectory where all the next steps will be. Summarise read quality and number of reads for each plate.

```
cd MiFish/

for K in {1..2}; do
  qiime demux summarize \
    --i-data demux_Plate$K.qza \
    --o-visualization demux_Plate$K.qzv
done

qiime demux summarize \
    --i-data demux_extras.qza \
    --o-visualization demux_extras.qzv

```

View ```.qzv``` files at [view.qiime2.org](view.qiime2.org). There is a dip in quality around 80-90bp. This could be due to primer dimers. Most samples have tens of thousands to hundreds of thousands of reads, which is more than enough.


## 2. Trim primers using cutadapt

The MiFish sequences are:   

F primer: GTCGGTAAAACTCGTGCCAGC (21 bp)   
R primer: CATAGTGGGGTATCTAATCCCAGTTTG (27 bp)   

### Trim 3' ends first
At the 3' end of the read, the primer will have been read through after reading the MiFish region. I need to be looking for the reverse complement of the reverse primer in read 1 (—p-adapter-f) and the reverse complement of the forward primer in R2 (—p-adapter-r).

F primer reverse complement: GCTGGCACGAGTTTTACCGAC   
R primer reverse complement: CAAACTGGGATTAGATACCCCACTATG    

```
for K in {1..2}; do
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux_Plate$K.qza \
    --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
    --p-adapter-r GCTGGCACGAGTTTTACCGAC \
    --o-trimmed-sequences trimd_Plate$K.qza \
    --verbose > cutadapt_out_Plate$K.txt
done

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux_extras.qza \
    --p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
    --p-adapter-r GCTGGCACGAGTTTTACCGAC \
    --o-trimmed-sequences trimd_extras.qza \
    --verbose > cutadapt_out_extras.txt
```

To see how much data passed the filter for each sample:
```
grep "Total written (filtered):" cutadapt_out_Plate1.txt 
grep "Total written (filtered):" cutadapt_out_Plate2.txt
grep "Total written (filtered):" cutadapt_out_extras.txt
```
Very variable amounts of data passed the filters, which is a shame.

Make new visualisations to see how many sequences are left and their quality scores.
```
for K in {1..2}; do
  qiime demux summarize \
    --i-data trimd_Plate$K.qza \
    --o-visualization trimd_Plate$K.qzv
done 

qiime demux summarize \
    --i-data trimd_extras.qza \
    --o-visualization trimd_extras.qzv
```
Looking at the trimd.qzv files, there are still a few thousand sequences for each sample at least.

### Trim 5' ends of reads
All R1 should begin with the forward primer: GTCGGTAAAACTCGTGCCAGC (21 bases).
All R2 should begin with the reverse primer: CATAGTGGGGTATCTAATCCCAGTTTG (27 bases).

Trim these with the following commands:

```
for K in {1..2}; do
  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences trimd_Plate$K.qza \
    --p-front-f GTCGGTAAAACTCGTGCCAGC \
    --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
    --o-trimmed-sequences trimd2_Plate$K.qza \
    --verbose > cutadapt_out2_Plate$K.txt
done

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences trimd_extras.qza \
  --p-front-f GTCGGTAAAACTCGTGCCAGC \
  --p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
  --o-trimmed-sequences trimd2_extras.qza \
  --verbose > cutadapt_out2_extras.txt
```  
About 90% of the data seems to pass this filter for both puffins and terns.

Not going to make the qzv files for now.

To see how much data passed the filter for each sample:
```
grep "Total written (filtered):" cutadapt_out2_Plate1.txt 
grep "Total written (filtered):" cutadapt_out2_Plate2.txt
grep "Total written (filtered):" cutadapt_out2_extras.txt
```
About 90% passed the filter here for all samples. That's pretty normal.

## 3. Denoise with dada2
I am going to use the same settings that I used for the 2017 and 2018 tern fecal samples here, except I need to add the --p-min-overlap parameter, otherwise I seem to be getting a load of rubbish reads which are 250bp long and start with long strings of Cs. This is only available in qiime2-2021.4 and later. I initally tried specifying an overlap of 30, but that didn't seem like enough as I was still getting junk sequences, but I think 50 is working well now.

Note, this step is pretty slow to run, a whole plate takes about 40 mins (but that changes depending on sequencing depth).
```
for K in {1..2}; do
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

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimd2_extras.qza \
  --p-trunc-len-f 133 \
  --p-trunc-len-r 138 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-min-overlap 50 \
  --p-n-threads 16 \
  --o-representative-sequences rep-seqs_extras \
  --o-table table_extras \
  --o-denoising-stats denoise_extras
```

Create visualizations for the denoising stats.
```
for K in {1..2}; do  
  qiime metadata tabulate\
    --m-input-file denoise_Plate$K.qza\
    --o-visualization denoise_Plate$K.qzv
done

qiime metadata tabulate\
  --m-input-file denoise_extras.qza\
  --o-visualization denoise_extras.qzv
```
Very variable numbers of reads get through the denoising/merging/chimera filtering. Not the best stats but hopefully it's fine.

To view the rep-seqs (only running this for one plate, will view all after merging plates)
```
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_Plate1.qza \
  --o-visualization rep-seqs_Plate1
```
This looks good. Using 50 as a minimum overlap seems to get rid of junk sequences.


## 4. Merge across plates
I need to merge both the feature tables, which contain the counts of each feature, and the rep-seqs, which contain the actual sequence for each feature. I am using the ```sum``` overlap method so that the data from the samples we sent in as tests are added to the data we got from them when we sequenced them again on the plate.

Note that this requires metadata for each sample, contained in metadata.txt. At the moment this is prety sparse as I didn't have time to add it all, but if you populate it with breeding stage, breeding period, substrate etc then it's easy to visalise those different groups in the later steps of this analysis, so I would recommend editing the metadata file. You can add any number of new columns to add info about the samples. It just has to be saved as a tab-separated file. 
```
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

```

## 5. Create a database for taxonomy assignment

The rescript module can download sequences from GenBank (NCBI) given an entrez search term. In this search, I'm searching for anything  with "12S", "12s", "small subunit ribosomal RNA gene", or whole mitochondrial genomes in the title. Then the rest of the terms should capture fish. However, I also want to include seabirds in the database, so that any sequences from the puffins get assigned, so I've added some terms for the birds I work with (so that I don't have to recreate this anew every time). 
I'm then filtering to just mitochondrial DNA in case any nuclear sequences get captured. 

```
qiime rescript get-ncbi-data \
--p-query '(12s[Title] OR \
            12S[Title] OR \
            "small subunit ribosomal RNA gene"[Title] OR \
            "mitochondrial DNA, complete genome"[Title] AND \
            ("Chondrichthyes"[Organism] OR \
            "Dipnomorpha"[Organism] OR \
            "Actinopterygii"[Organism] OR \
            "Myxini"[Organism] OR \
            "Hyperoartia"[Organism] OR \
            "Coelacanthimorpha"[Organism] OR \
            fish[All Fields] OR \
            "Sphenisciformes"[Organism] OR \
            "Charadriiformes"[Organism] OR \
            "Procellariiformes[Organism])) \
            AND mitochondrion[filter]' \
  --p-n-jobs 5 \
  --o-sequences ncbi-refseqs-unfiltered.qza \
  --o-taxonomy ncbi-taxonomy-unfiltered.qza  
```

The next step is to clean up the database. They recommend removing sequences with too many degenerate bases or long runs of homoploymers using the cull-seqs option. I am running with default params, which will remove sequences with more than 5 degenerate bases and runs of 8 homopolymers.
```
qiime rescript cull-seqs \
  --i-sequences ncbi-refseqs-unfiltered.qza \
  --p-num-degenerates 5 \
  --p-homopolymer-length 8 \
  --p-n-jobs 4 \
  --o-clean-sequences ncbi-refseqs-culled
```

The final step is to dereplicate the database.
```
qiime rescript dereplicate \
  --i-sequences ncbi-refseqs-culled.qza \
  --i-taxa ncbi-taxonomy-unfiltered.qza \
  --p-mode uniq \
  --p-threads 4 \
  --o-dereplicated-sequences ncbi-refseqs-culled-derep \
  --o-dereplicated-taxa ncbi-taxonomy-culled-derep
```

I think this step also takes care of cleaning up the taxonomy file, so if any sequences were culled in the above step, then they should be removed.

In the summer I was playing around with training a feature classifier on a similar database to assign taxonomy, but it did really badly. Butterfish were consistently assigned to the wrong species, as were haddock I think. So I will use the script I have to blast all the rep-seqs against the reference database.

Note the blast method only works with an older version of Qiime, so I have to load that first. 

## 6. Assign taxonomy

In the summer I was playing around with training a feature classifier on a similar database to assign taxonomy, but it did really badly. Butterfish were consistently assigned to the wrong species, as were haddock I think. So I will use the script I have to blast all the rep-seqs against the reference database.

Note the blast method only works with an older version of Qiime, so I have to load that first. 
```
conda activate qiime2-2019.4

./mktaxa.py ncbi-refseqs-culled-derep.qza ncbi-taxonomy-culled-derep.qza merged_rep-seqs.qza
```
This creates a file called ```superblast_taxonomy.qza```, which relates the rep-seqs to species. Make a version that can be viewed online:
```
qiime metadata tabulate \
  --m-input-file superblast_taxonomy.qza \
  --o-visualization ncbi-refseqs-superblast_taxonomy
```
This file can be used to look up the assignments made to individual features, and to search for species e.g. looking up "Aves" shows that nearly all sequences assigned to Aves were classified as puffins. 

### Add a human sequence to the database manually
I made an error and forgot to add humans to the database. I don't want to add them to my search term, as it will download thousands of human sequences, so I just found one complete mitochondrial genome on GenBank. I put this into a text document ```temp.txt``` and deleted the carriage returns at the end of each line using ```tr -d '\n' < temp.txt > output.txt ```. 

I exported the sequences and taxonomy artefacts from qiime into simple text docs:
```
qiime tools export \
  --input-path ncbi-taxonomy-culled-derep.qza \
  --output-path editing_database
  
qiime tools export \
  --input-path ncbi-refseqs-culled-derep.qza \
  --output-path editing_database
  
```
I added the sequence and corresponding taxonomy string to each file. Then I imported them back into Qiime:
```
qiime tools import \
  --input-path editing_database/taxonomy.tsv \
  --output-path editing_database/ncbi-taxonomy-withHuman \ 
  --type 'FeatureData[Taxonomy]'
  
qiime tools import \
  --input-path editing_database/dna-sequences.fasta \
  --output-path editing_database/ncbi-refseqs-withHuman \
  --type 'FeatureData[Sequence]'  
```
### Re-assign taxonomy using edited database
Note, I am just overwriting the old file here since I don't want to have two versions hanging around.
```
conda activate qiime2-2019.4

./mktaxa.py editing_database/ncbi-refseqs-withHuman.qza editing_database/ncbi-taxonomy-withHuman.qza merged_rep-seqs.qza
```
Remake the version that can be viewed online:
```
qiime metadata tabulate \
  --m-input-file superblast_taxonomy.qza \
  --o-visualization ncbi-refseqs-superblast_taxonomy
```
Searching this shows five features that were classified as human, although three of them have <95% identity, so are likely some other mammal. 

## 7. Remove non-food reads
Filter out any sequences from the bird, mammals (human), and unnassigned sequences before rarefying. 

```
qiime taxa filter-table \
  --i-table merged-table.qza \
  --i-taxonomy superblast_taxonomy.qza \
  --p-exclude Unassigned,Aves,Mammalia \
  --o-filtered-table merged_table_noBirdsMammalsUnassigned.qza
  
qiime feature-table summarize \
    --i-table merged_table_noBirdsMammalsUnassigned.qza \
    --m-sample-metadata-file metadata.txt \
    --o-visualization merged_table_noBirdsMammalsUnassigned
```

## 8. Make some barplots
First, I do actually want to see how many puffin, human, and unassigned sequences there were.
```
qiime taxa barplot \
  --i-table merged-table.qza \
  --i-taxonomy superblast_taxonomy.qza \
  --m-metadata-file metadata.txt \
  --o-visualization barplot_before_filtering.qzv
```
19 samples had 100% avian DNA, which is annoying since we can't recover any fish DNA from them. I haven't seen this in terns. Perhaps the primers are a better match to puffins. Or perhaps the birds hadn't fed on anything recently?

Barplot having removed bird/human/unassigned DNA:
```
qiime taxa barplot \
  --i-table merged_table_noBirdsMammalsUnassigned.qza \
  --i-taxonomy superblast_taxonomy.qza \
  --m-metadata-file metadata.txt \
  --o-visualization barplot_noBirdsMammalsUnassigned.qzv
```

## 9. Calculate alpha rarefaction curves 
First you have to collapse the taxonomy, so that you are rarefying based on species assignments and not ASVs. I am going to rarefy from a depth of 100 to 20,000 reads.
```
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
```
Select sample metadata column: SampleType to separate the samples out from the blanks and mock community. The Shannon diversity is flat for the samples (ignore the mock community, this is artifically diverse compared to the samples). The curve for the observed OTUs (species) increases from 100 to 2000, so repeat the rarefaction for this range to find out what the actiual minimum depth is without sacrificing the diversity of the samples.
```
qiime diversity alpha-rarefaction \
  --i-table merged_table_noBirdsMammalsUnassigned_collapsed.qza \
  --m-metadata-file metadata.txt \
  --p-min-depth 100 \
  --p-max-depth 2500 \
  --o-visualization alpha-rarefaction-100-2500
```
Shannon diversity is still flat for the samples. The observed OTUs has quite a lot of variability at 400 but has plateued by 600. 