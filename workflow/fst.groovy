/*===========================================================

	P I P E L I N E   P A R A M E T E R S

===========================================================*/

fst_program = new File ("..").getCanonicalPath() // location of this program, i.e. /path/to/program/, default is level up from this script
fst_output_folder = "$fst_program/tests/SFPQ-ABL1/output" // where all the results will live

/*--Inputs--*/
fusion_finder_output = "$fst_program/tests/SFPQ-ABL1/input/jaffa_results.csv" // Location of fusion finders output file
column_positions = "3 4 5 6" // Location of chromosome 1, break point 1, chromosome 2, breakpoint 2, columns in above file
delimiter = "c" // is the above file tab or comma delimited (t or c, c if ommited)
genome = "19" // was the above file generated with hg19 or hg38 (19, 38)

/*--Alignment Parameters--*/
fastq_file_directory = "/mnt/CANC1-genomic/RNA_Seq_Experiments/Fastq_files/20160610_20160722/"
threads = "16" // How many threads should the aligner use

/*--Print a fusion of interest?--*/
find_fusion = true //true will create generate a plot of below fusion, false will stop after the alignment.
fusion = "SFPQ:ABL1" //if true above, identify a fusion of interest. Must be in order of fusion and seperated with a colon (:)"
fusion_friendly = fusion.replaceAll(":", "_") //allows for windows compatiability

/*=======================================================

	P I P E L I N E   P A T H S
	Note: change only if comfortable editing R (fst_plot.r)

=======================================================*/

annotation_folder = "$fst_output_folder/annotation"
alignment_folder = "$fst_output_folder/alignment/pass1"
genome_folder = "$fst_output_folder/genome"
reference_folder = "$fst_output_folder/reference"
fusion_folder = "$alignment_folder/$fusion_friendly"

/*=======================================================

	P I P E L I N E   S T A G E S

=======================================================*/

/*--Generate fused annotation and fst reference fasta files--*/

generate_fst = {
	exec "python $fst_program/fusion/main.py bpipe -in $fusion_finder_output -out $fst_output_folder -pos $column_positions -del $delimiter -genome $genome"
	exec "python $fst_program/fusion/exon_annotation.py $annotation_folder/jaffa_annotation.bed $annotation_folder/exon_boundaries.bed $annotation_folder/gene_boundaries.bed"
}

/*--Align to the fst reference fasta created in previous stage
	NOTE: Stop after this stage to view in IGV--*/

star_align = {
	exec "STAR --runMode genomeGenerate --genomeDir $genome_folder --genomeFastaFiles $reference_folder/fst_reference.fasta --runThreadN $threads --limitGenomeGenerateRAM 32000000000"
	exec "STAR --genomeDir $genome_folder --readFilesIn $fastq_file_1 $fastq_file_2 --runThreadN $threads --outWigType wiggle --outWigNorm RPM --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix $alignment_folder/ --limitBAMsortRAM 10000000000"
}

/*--Prepare fusion files for visualisation--*/

prepare_fusion_plot = {
	exec "$fst_program/plotit/fst_plot_prep.sh $fusion $fusion_folder $fusion_friendly $annotation_folder $alignment_folder $reference_folder"
}

/*--Print fusion pdf--*/

plot_fusion = {
	exec "Rscript $fst_program/plotit/fst_plot.r $fst_output_folder $fusion"
}

/*--Run bpipe, run--*/

if(find_fusion){
	Bpipe.run {
		generate_fst + prepare_fusion_plot + plot_fusion
		//generate_fst + star_align + prepare_fusion_plot + plot_fusion
	}
} else {
	Bpipe.run {
		generate_fst + star_align
	}
}
