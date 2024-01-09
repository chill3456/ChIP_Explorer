#!/bin/bash 

#applications needed 
#need to install bedtools https://bedtools.readthedocs.io/en/latest/content/installation.html

#need to install homer http://homer.ucsd.edu/homer/introduction/install.html 

#need to install MEME-ChIP https://meme-suite.org/meme/doc/install.html?man_type=web


#files needed 

#TSS cordinates with gene name in second to last column like refGene_hg19_TSS_srt.bed can be downloaded from here https://sourceforge.net/projects/seqminer/files/Reference%20coordinate/refGene_hg19_TSS.bed/download

#should look like this 

#chr1	11873	11873	NR_046018	DDX11L1	+

#hg19.fa or you genomes fasta file can been downloaded from ucsc https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/






# Initialize variables
output_folder=""
declare -a peak_files

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --peaks)
            shift # past the keyword
            while [[ $# -gt 0 ]] && [[ ! $1 == *'--'* ]]; do
                peak_files+=("$1") # add to array
                shift # past the value
            done
            ;;
        --output_folder)
            output_folder="$2"
            shift # past argument
            shift # past value
            ;;
        *) # unknown option
            shift # past argument
            ;;
    esac
done

# Check if output folder is specified
if [ -z "$output_folder" ]; then
    echo "Output folder not specified."
    exit 1
fi

# Create the Peak_analysis subfolder
peak_analysis_folder="$output_folder/Peak_analysis/"
mkdir -p "$peak_analysis_folder"

# Process each file
for file in "${peak_files[@]}"; do
    # Make sure file exists
    if [ -f "$file" ]; then
        echo "Sorting $file..."

        # Construct the output file name
        base_name=$(basename "$file" .bed)
        sorted_file="$peak_analysis_folder/${base_name}_sorted.bed"

        # Run bedtools sort and save the output
        bedtools sort -i "$file" > "$sorted_file"

    else
        echo "File $file not found."
    fi
done

cd $output_folder/Peak_analysis/

mkdir enhancers_and_promoters_peaks

for i in *_sorted.bed
	do 
		#assign names of files have to rename for out files 
		NARROWPEAK_IN=$i

		NEAREST_GENE=${i/_sorted.bed/_nearest_gene.bed}

			bedtools closest -a $NARROWPEAK_IN -b /pathto/refGene_hg19_TSS_sorted.bed -d -t first > $output_folder/Peak_analysis/enhancers_and_promoters_peaks/$NEAREST_GENE

		done

cd $output_folder/Peak_analysis/enhancers_and_promoters_peaks





for k in *_nearest_gene.bed 
	do 
		#assign names of files have to rename for out files 
		NEAREST_GENE_file=$k
		Enhancers_out=${k/_nearest_gene.bed/_enhancers_nearest_gene.bed}
		Promoters_out=${k/_nearest_gene.bed/_promoters_nearest_gene.bed}
		All_out=${k/_nearest_gene.bed/_all_nearest_gene.bed}
			#finding enhancers (greater then 1000bp noted in last column) and promoters
			awk ' $(NF) > 1000 ' $NEAREST_GENE_file > $output_folder/Peak_analysis/enhancers_and_promoters_peaks/$Enhancers_out
			awk ' $(NF) <= 1000 ' $NEAREST_GENE_file > $output_folder/Peak_analysis/enhancers_and_promoters_peaks/$Promoters_out
			mv $NEAREST_GENE_file $All_out

		done

#MAKE TABLE OF PEAK COUNTS FOR ALL SAMPLES
for i in *nearest_gene.bed; do wc -l $i; done | sed s/_nearest_gene.bed//g | awk -v OFS="\t" '{print $2,$1}' > $output_folder/Peak_analysis/enhancers_and_promoters_peaks/nearest_gene_counts.txt

cd $output_folder/Peak_analysis/

mkdir peak_overlaps_enhancers 

mkdir peak_overlaps_enhancers_not

mkdir peak_overlaps_all

mkdir peak_overlaps_all_not

mkdir peak_overlaps_promoters

mkdir peak_overlaps_promoters_not



cd $output_folder/Peak_analysis/enhancers_and_promoters_peaks

#doing peak laps for all peaks 
for FIRST in $(find *_all_nearest_gene.bed)
do
for SECOND in $(find *_all_nearest_gene.bed*)
do
		#do pairwise comparision first finding non overlapping peaks
		bedtools intersect -a $FIRST -b $SECOND -v > $output_folder/Peak_analysis/peak_overlaps_all_not/${FIRST}_NOT_overlaping_${SECOND}_overlaping.bed
		#then overlapping peaks
		bedtools intersect -a $FIRST -b $SECOND > $output_folder/Peak_analysis/peak_overlaps_all/${FIRST}_overlaping_${SECOND}_overlaping.bed
done
done


#doing peak laps for promoters 
for FIRST in $(find *_promoters_nearest_gene.bed)
do
for SECOND in $(find *_promoters_nearest_gene.bed)
do
		#do pairwise comparision ffirst finding non overlapping peaks
		bedtools intersect -a $FIRST -b $SECOND -v > $output_folder/Peak_analysis/peak_overlaps_promoters_not/${FIRST}_NOT_overlaping_${SECOND}_overlaping.bed
		#then overlapping peaks
		bedtools intersect -a $FIRST -b $SECOND > $output_folder/Peak_analysis/peak_overlaps_promoters/${FIRST}_overlaping_${SECOND}_overlaping.bed
done
done


#doing peak laps for enhancers
for FIRST in $(find *_enhancers_nearest_gene.bed)
do
for SECOND in $(find *_enhancers_nearest_gene.bed)
do
		#do pairwise comparision ffirst finding non overlapping peaks
		bedtools intersect -a $FIRST -b $SECOND -v > $output_folder/Peak_analysis/peak_overlaps_enhancers_not/${FIRST}_NOT_${SECOND}_overlaping.bed
		#then overlapping peaks
		bedtools intersect -a $FIRST -b $SECOND > $output_folder/Peak_analysis/peak_overlaps_enhancers/${FIRST}_overlaping_${SECOND}_overlaping.bed
done
done


cd $output_folder/Peak_analysis/peak_overlaps_enhancers_not


#MAKE TABLE OF PEAK COUNTS FOR ALL SAMPLES
for i in *overlaping.bed; do wc -l $i; done | sed s/_overlaping.bed//g | awk -v OFS="\t" '{print $2,$1}' > $output_folder/Peak_analysis/peak_overlaps_enhancers_not/overlaping_peak_counts.txt

cd $output_folder/Peak_analysis/peak_overlaps_enhancers

#MAKE TABLE OF PEAK COUNTS FOR ALL SAMPLES
for i in *overlaping.bed; do wc -l $i; done | sed s/_overlaping.bed//g | awk -v OFS="\t" '{print $2,$1}' > $output_folder/Peak_analysis/peak_overlaps_enhancers/overlaping_peak_counts.txt

cd $output_folder/Peak_analysis/peak_overlaps_promoters_not

#MAKE TABLE OF PEAK COUNTS FOR ALL SAMPLES
for i in *overlaping.bed; do wc -l $i; done | sed s/_overlaping.bed//g | awk -v OFS="\t" '{print $2,$1}' > $output_folder/Peak_analysis/peak_overlaps_promoters_not/overlaping_peak_counts.txt

cd $output_folder/Peak_analysis/peak_overlaps_promoters

#MAKE TABLE OF PEAK COUNTS FOR ALL SAMPLES
for i in *overlaping.bed; do wc -l $i; done | sed s/_overlaping.bed//g | awk -v OFS="\t" '{print $2,$1}' > $output_folder/Peak_analysis/peak_overlaps_promoters/overlaping_peak_counts.txt

cd $output_folder/Peak_analysis/peak_overlaps_all_not

#MAKE TABLE OF PEAK COUNTS FOR ALL SAMPLES
for i in *overlaping.bed; do wc -l $i; done | sed s/_overlaping.bed//g | awk -v OFS="\t" '{print $2,$1}' > $output_folder/Peak_analysis/peak_overlaps_all_not/overlaping_peak_counts.txt

cd $output_folder/Peak_analysis/peak_overlaps_all

#MAKE TABLE OF PEAK COUNTS FOR ALL SAMPLES
for i in *overlaping.bed; do wc -l $i; done | sed s/_overlaping.bed//g | awk -v OFS="\t" '{print $2,$1}' > $output_folder/Peak_analysis/peak_overlaps_all/overlaping_peak_counts.txt


#now doing motif analysis of enhancers_and_promoters_peaks



#HOMER
cd $output_folder/Peak_analysis/enhancers_and_promoters_peaks

PATH=$PATH:/pathto/homer/bin


for l in *.bed; do  Overlaping10_file=$l HOMER_out=${l/.bed/_HOMER}; findMotifsGenome.pl $Overlaping10_file hg19 $HOMER_out -preparse -size 200 -mask -preparsedDir /pathto/homer/data/genomes/hg19; done
		

#MEME-ChIP




for l in *.bed
  do 
    #assign names of files have to rename for out files 
    bed=$l
    FASTA_out=${l/.bed/.fasta}
      #adding 
      bedtools getfasta -fi /pathto/hg19.fa -bed $bed -fo $FASTA_out

    done




for l in *.fasta; do  FASTA=$l MEME_out=${l/.fasta/_MEME_ChIP}; meme-chip -oc . -time 240 -ccut 100 -dna -order 2 -minw 6 -maxw 15 -db /pathto/meme/db/motif_databases/MOUSE/uniprobe_mouse.meme -db /pathto/meme/db/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme -db /pathto/meme/db/motif_databases/EUKARYOTE/jolma2013.meme -meme-mod zoops -meme-nmotifs 3 -meme-searchsize 100000 -streme-pvt 0.05 -streme-totallength 4000000 -centrimo-score 5.0 -centrimo-ethresh 10.0 $FASTA -oc $MEME_out; done



echo "Processing completed. Output files are in $output_folder."