#!/bin/bash
#------------------------------------------------------------------------------------------
# Get variables 
#------------------------------------------------------------------------------------------
current_directory=`dirname $0`
fusion=$1
fusion_folder=$2
fusion_friendly=$3
annotation_folder=$4
alignment_folder=$5
reference_folder=$6
fst_program=$7
#------------------------------------------------------------------------------------------
# Make a fusion folder directory
#------------------------------------------------------------------------------------------
if [ ! -d $fusion_folder ];
	then mkdir $fusion_folder;
fi
#------------------------------------------------------------------------------------------
# Reduce BAM file to include only reads from desired fusion, this will reduce space
#------------------------------------------------------------------------------------------
printf '\n'
echo $fusion
echo '------------------------------------------'
echo 'filtering BAM file for fusion of interest'
samtools view -hq 1 $alignment_folder/Aligned.sortedByCoord.out.bam {$fusion} | \
			samtools view -Sb - > $fusion_folder/${fusion_friendly}.bam
samtools index $fusion_folder/${fusion_friendly}.bam
#------------------------------------------------------------------------------------------
# Filter out reads with overhangs less than 5 (accounts for much of the noise in final plot)
#------------------------------------------------------------------------------------------
echo 'filtering BAM file for reads with overhangs < 2 (noise reduction)'
samtools view -h $fusion_folder/${fusion_friendly}.bam | \
awk '{if($0 ~ /^@/){print $0;} else {cigar = $6; gsub("([0-9]+[S])","",cigar); gsub("([0-9]+[I])","",cigar); gsub("([0-9]+[D])","",cigar); gsub("([0-9]+[N])"S,"",cigar); m = split(cigar,matches,"M"); if(m > 2) {if(matches[1] > 0 && matches[m-1] > 0){print $0;}} else if(m <= 2) {print $0;}}}' | \
samtools view -Sb - > $fusion_folder/${fusion_friendly}_lt5.bam
#------------------------------------------------------------------------------------------
# Create files including only reads involved in fusion and only the other reads
#------------------------------------------------------------------------------------------
echo 'Creating ancillilary files'
grep $fusion $annotation_folder/gene_boundaries.bed | awk '$5 == 1' | sort | head -n 1 | awk -v fusion="$fusion" '{print fusion"\t"$2"\t"$2+1}' > $fusion_folder/fusion_locations.bed
samtools view $fusion_folder/${fusion_friendly}_lt5.bam -b -h -o $fusion_folder/fusion_reads.bam -U $fusion_folder/splice_junction_reads.bam -L $fusion_folder/fusion_locations.bed
samtools view -h $fusion_folder/fusion_reads.bam | awk '{if($0 ~ /^@/ || $6 ~ /N/) {print $0}}' | samtools view -Sb - > $fusion_folder/fusion_breakpoint_reads.bam
#------------------------------------------------------------------------------------------
# Index the new bam files
#------------------------------------------------------------------------------------------
echo 'Index BAM files'
samtools index $fusion_folder/splice_junction_reads.bam
samtools index $fusion_folder/fusion_breakpoint_reads.bam
#------------------------------------------------------------------------------------------
# Remove Junk
#------------------------------------------------------------------------------------------
rm $fusion_folder/fusion_reads.bam
#------------------------------------------------------------------------------------------
# Create SJ.Tab.out
#------------------------------------------------------------------------------------------
samtools view -h -o $fusion_folder/${fusion_friendly}.sam $fusion_folder/${fusion_friendly}_lt5.bam
awk -f $current_directory/sj_out_gen.awk $fusion_folder/${fusion_friendly}.sam | sort > $fusion_folder/junctions.txt
