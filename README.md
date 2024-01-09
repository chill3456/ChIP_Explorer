<p align="center">
  <a>
    <img width="40%" src="https://github.com/chill3456/ChIP_Explorer/blob/master/assets/graphic.png" alt="graphic describing ChIP_Explorer.sh png">
  </a>
  <h1 align="center">ChIP_Explorer</h1>
</p>




# Description 

ChIP_Exploration.R explores of ChIP-seq results by annotating peaks as enhancers or promoters, finding peak overlaps, and motif analysis of peaks in single Bash script that can be run on the command line. The script is configured for the hg19 assembly of thehuman genome but can easily be adjusted for the required genome.

# Required applications 

`bedtools`
`homer`
`MEME-ChIP`

# Required files

Files for hg19 are provided and can also be downloaded from https://sourceforge.net/projects/seqminer/files/Reference%20coordinate/refGene_hg19_TSS.bed/download and https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/. They are a bed file containing all the transcription start sites (TSS) in the genome and a fasta file for the hg19 genome. If not using hg19 other equivalent files for other genomes can be downloaded.

`required_files/refGene_hg19_TSS_sorted.bed`

`required_files/hg19.fa`

# Input files 
Bed files with the first three columns being: chromosome, PeakChromosomeStartPosition, PeakChromosomeEndPosition. First row from example bed file /example_files/peaksSet3.bed

`chr10	69644114	69644598`

# Running the script 

The script is a BASH script and can be run on the command line using sbatch or sh. 

`sbatch ChIP_Explorer.sh --peaks /example_files/peaksSet1.bed /example_files/peaksSet2.bed /example_files/peaksSet3.bed --output_folder /output_examples/`

# Enhancer Promoter Outputs 

The script annotates each peak in the .bed file as either an enhancer or promoter site based on its distance to the nearest transcription start site (TSS) in the refGene_hg19_TSS_sorted.bed. Peaks are annotated as enhancers if they are more than 1kb from the nearest TSS and promoters if they are 1kb or less from the nearest TSS. Bed files containing all peaks, just enhancers, and just promoters are output. 

`/outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_all_nearest_gene.bed`

`/outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_enhancers_nearest_gene.bed`

`/outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_promoters_nearest_gene.bed`

These files contain the gene name of the gene which TSS the peaks was closest to in the third to last column, the strand of this gene in the second to last column, and the distance from the peak to nearest tss in the 3rd to last column. First row from example bed file /outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_enhancers_nearest_gene.bed. 

`chr1	10016	10491	chr1	11873	11873	NR_046018	DDX11L1	+	1382`

A summary file that displays the amount of peaks in each all peaks, enhancer, and promoter file is output /outputExamples/Peak_analysis/enhancers_and_promoters_peaks/nearest_gene_counts.txt. 

| peaksSet1_all	      | 4051 | 
| peaksSet1_enhancers	| 1946 | 
| peaksSet1_promoters	| 2105 | 
| peaksSet2_all       | 226  | 
| peaksSet2_enhancers	| 181  | 
| peaksSet2_promoters	| 45   | 
| peaksSet3_all	      | 30   | 
| peaksSet3_enhancers	| 13   | 
| peaksSet3_promoters	| 17   | 

# Peak Overlap Outputs 

The script finds overlap of peaks between the different data sets and generates bed files contain peaks that overlap or dont overlap. Pairwise comparisons are done between all the peaks, just enhancers, or just promoters. Folders contain each type of pairwise comparison with either overlaping of not overlapping peaks (those ending in _not). Folders generated:

`peak_overlaps_all`

`peak_overlaps_all_not`

`peak_overlaps_enhancers`

`peak_overlaps_enhancers_not`

`peak_overlaps_promoters`

`peak_overlaps_promoters_not`

Bed files are named based on the comparison with first peak set followed by _overlaping or _NOT_overlaping and then the peak set it was examine for overlap or not overlap with. They only contain peaks from the first peak set. 

`peaksSet1_enhancers_nearest_gene.bed_overlaping_peaksSet2_enhancers_nearest_gene.bed_overlaping.bed`

`peaksSet1_promoters_nearest_gene.bed_NOT_overlaping_peaksSet2_promoters_nearest_gene.bed_overlaping.bed`

A summary file that displays the amount of peaks in each type of comparisons is output /outputExamples/Peak_analysis/peak_overlaps_all/overlaping_peak_counts.txt


| peaksSet1_all_nearest_gene.bed_overlaping_peaksSet1_all_nearest_gene.bed | 4051 | 
| peaksSet1_all_nearest_gene.bed_overlaping_peaksSet2_all_nearest_gene.bed | 22   | 
| peaksSet1_all_nearest_gene.bed_overlaping_peaksSet3_all_nearest_gene.bed | 18   | 
| peaksSet2_all_nearest_gene.bed_overlaping_peaksSet1_all_nearest_gene.bed | 22   | 
| peaksSet2_all_nearest_gene.bed_overlaping_peaksSet2_all_nearest_gene.bed | 226  | 
| peaksSet2_all_nearest_gene.bed_overlaping_peaksSet3_all_nearest_gene.bed | 1    | 
| peaksSet3_all_nearest_gene.bed_overlaping_peaksSet1_all_nearest_gene.bed | 18   | 
| peaksSet3_all_nearest_gene.bed_overlaping_peaksSet2_all_nearest_gene.bed | 1    | 
| peaksSet3_all_nearest_gene.bed_overlaping_peaksSet3_all_nearest_gene.bed | 30   | 


# Motif Analysis Output

The final output is motif analysis is done on all BED files in the /outputExamples/Peak_analysis/enhancers_and_promoters_peaks folder containing all peaks, just enhancers, and just promoters. Motif analysis is done using both HOMER and MEME-ChIP. 

HOMER results are output in a folder with the BED file name and ending in _HOMER such as /outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_all_nearest_gene_HOMER. 

The output contains /outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_all_nearest_gene_HOMER/knownResults.html that looks for motifs found in previous ChIP-seqs.

<img src='https://github.com/chill3456/Quick_qPCR_Analysis/blob/master/assets/knownResults.png' width = '100%' alt="knownResults png">

The output also contains /outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_all_nearest_gene_HOMER/homerResults.html that looks for de-novo motifs and then looks for similarity to established motifs.

<img src='https://github.com/chill3456/Quick_qPCR_Analysis/blob/master/assets/homerResults.png' width = '100%' alt="homerResults png">

MEME-ChIP results are output in a folder with the BED file name and ending in _MEME_ChIP such as /outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_all_nearest_gene_MEME_ChIP. 

The output is summarized in /outputExamples/Peak_analysis/enhancers_and_promoters_peaks/peaksSet1_all_nearest_gene_MEME_ChIP/meme-chip.html. 

<img src='https://github.com/chill3456/Quick_qPCR_Analysis/blob/master/assets/meme-chip.png' width = '100%' alt="meme-chip png">

