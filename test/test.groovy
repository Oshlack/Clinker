/*===========================================================

	P I P E L I N E   P A R A M E T E R S

===========================================================*/

fst_program = new File ("..").getCanonicalPath() // location of this program, i.e. /path/to/program/, default is level up from this script
fst_output_folder = "../results/test" // where all the results will live

/*--Inputs--*/
fusion_finder_output = "caller/bcr_abl1.csv" // Location of fusion finders output file
column_positions = "1 2 3 4" // Location of chromosome 1, break point 1, chromosome 2, breakpoint 2, columns in above file
delimiter = "c" // is the above file tab or comma delimited (t or c, c if ommited)
genome = "19" // was the above file generated with hg19 or hg38 (19, 38)
competitive = true

/*--Alignment Parameters--*/
threads = "16" // How many threads should the aligner use

/*--Print a fusion of interest?--*/
find_fusion = true //true will create generate a plot of below fusion, false will stop after the alignment.
fusion = "BCR:ABL1" //if true above, identify a fusion of interest. Must be in order of fusion and seperated with a colon (:)"
fusion_friendly = fusion.replaceAll(":", "_") //allows for windows compatiability

pdf_width = "8.27"
pdf_height = "11.69"
proportion = "1,3,1,3,3,2"

/*=======================================================

	P I P E L I N E   P A T H S
	Note: change only if comfortable editing R (fst_plot.r)

=======================================================*/

annotation_folder = "$fst_output_folder/annotation"
alignment_folder = "$fst_output_folder/alignment"
genome_folder = "$fst_output_folder/genome"
reference_folder = "$fst_output_folder/reference"
fusion_folder = "$alignment_folder/$fusion_friendly"

/*=======================================================

	P I P E L I N E   S T A G E S

=======================================================*/

/*--Generate fused annotation and fst reference fasta files--*/

generate_fst = {
	produce("$fst_output_folder/reference/fst_reference.fasta") {
		exec "python $fst_program/fusion/main.py -in $fusion_finder_output -out $fst_output_folder -pos $column_positions -del $delimiter -genome $genome -competitive $competitive"
	}
}

star_genome_gen = {
    //Generate STAR genome index
    doc "Generate STAR genome index"

    produce("$genome_folder/Genome") {
	    exec """
	        STAR --runMode genomeGenerate
				--runThreadN $threads
				--genomeDir $genome_folder
				--genomeFastaFiles $reference_folder/fst_reference.fasta
				--limitGenomeGenerateRAM 32000000000
				--genomeSAindexNbases 5""","stargenind"
    }
}

star_align = {

	doc "Map paired-end reads using the STAR aligner: 1st pass"

  //Map paired-end reads using the STAR aligner: 1st pass
	from("*.fastq.gz") {
		transform("(.*)_R1.fastq.gz","(.*)_R2.fastq.gz") to ("Aligned.sortedByCoord.out.bam") {

			String R1 = inputs[0];
			String R2 = inputs[1];
			String files = R1+" "+R2;

			String file_name = R1.split("/")[-1]
			String folder_name = file_name.split("_R1")[0]
			new File("$alignment_folder/$folder_name").mkdirs()
			output.dir = "$alignment_folder/$folder_name/"

			exec """STAR --genomeDir $genome_folder
				--readFilesIn ${files}
				--readFilesCommand zcat
				--outSAMtype BAM SortedByCoordinate
				--runThreadN $threads
				--outFileNamePrefix "$output.dir/"
				--limitBAMsortRAM 20000000000
			""","star1pass"
	    }
	}
}

/*--Prepare fusion files for visualisation--*/

prepare_fusion_plot = {
	from("*.fastq.gz") {
		transform("(.*)_R1.fastq.gz","(.*)_R2.fastq.gz") to ("coverage_rpm.bedgraph") {

			String R1 = inputs[0];
			String file_name = R1.split("/")[-1]
			String folder_name = file_name.split("_R1")[0]
			output.dir = "$alignment_folder/$folder_name/$fusion_friendly"
			String new_alignment = "$alignment_folder/$folder_name"

			exec "$fst_program/plotit/fst_plot_prep.sh $fusion $output.dir $fusion_friendly $annotation_folder $new_alignment $reference_folder"
		}
	}
}

/*--Print fusion pdf--*/

plot_fusion = {
	from("*.fastq.gz") {
			String R1 = inputs[0];
			String file_name = R1.split("/")[-1]
			String folder_name = file_name.split("_R1")[0]
			output.dir = "$alignment_folder/$folder_name/$fusion_friendly"

			exec "Rscript $fst_program/plotit/fst_plot.r $fst_output_folder $fusion $output.dir $pdf_height $pdf_width $proportion"
	}
}

/*--Run bpipe, run--*/

if(find_fusion){
	Bpipe.run {
		generate_fst + star_genome_gen + "%_*.fastq.gz"*[star_align + prepare_fusion_plot + plot_fusion]
	}
} else {
	Bpipe.run {
		generate_fst + star_genome_gen + "%_*.fastq.gz"*[star_align]
	}
}
