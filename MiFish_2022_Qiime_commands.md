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

To see how much data passed the filter for each sample:

    grep "Total written (filtered):" cutadapt_out_Plate1.txt 
    grep "Total written (filtered):" cutadapt_out_Plate2.txt
    grep "Total written (filtered):" cutadapt_out_Plate3.txt
