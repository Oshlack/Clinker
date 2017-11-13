#===========================================================
#
# 	FST PRINTER
#
# 	Authors: Breon Schmidt & Anthony Hawkins
# 	Date: 29 Nov 2016
#
# 	NOTE: Excuse my R...
#
#===========================================================

#===========================================================
#
#   L I B R A R I E S
#
#===========================================================

suppressMessages(library(Gviz))
suppressMessages(library(IRanges))
suppressMessages(library(data.table))
options(ucscChromosomeNames=FALSE)

#===========================================================
#
#   F U N C T I O N S
#
#===========================================================


#---------------------------------------------------
#
#   Colour Generator
#
#   Input: Number of randomised colours required
#   Output: Vector of hex codes
#
#---------------------------------------------------

randomColour <- function(number){

	# Create storage to hold the colours
	colour_storage <- c(1:number)

	# Generate Hex codes, assign to storage
	for (i in 1:number)
	{
		hex_code <- sample(c(0:9, letters[1:6]), 6, replace=TRUE)
		colour_storage[i] <- paste(hex_code, collapse="", sep="")
		colour_storage[i] <- paste("#", colour_storage[i], sep="")
	}

	# Return list of colours
	return(colour_storage)
}

#---------------------------------------------------
#
#   Handle User Input Arguments,
#   including managing error checks
#
#   Input: N/A
#   Output: Location of results folder, fusion
#
#---------------------------------------------------

userInput <- function(){

	# Read in command line
	command_line <- commandArgs(trailingOnly = TRUE)

	# Unpack arguments
	results_location <- command_line[1]
	fusion <- as.character(command_line[2])
	fusion_location <- command_line[3]
	pdf_width <- command_line[4]
	pdf_height <- command_line[5]
	ratio <- command_line[6]
	alignment_folder <- command_line[7]
	sample_name <- command_line[8]
	track_colours <- command_line[9]
	track_order <- command_line[10]

	#results_location <- "Path to results folder"
	#fusion <- as.character("Your fusion name")

	# Return multiple variables
	multiple_return <- list(
		results_location = results_location, 
		fusion = fusion, 
		fusion_location = fusion_location, 
		pdf_width = pdf_width, 
		pdf_height = pdf_height, 
		ratio = ratio, 
		alignment_folder = alignment_folder, 
		sample_name = sample_name, 
		track_colours = track_colours,
		track_order = track_order
	)

	return(multiple_return)
}

#---------------------------------------------------
#
#   Track Width
#
#   Input: N/A
#   Output: Start and End of alignment/annotation tracks
#
#---------------------------------------------------

trackDimensions <- function(gene_file, fusion){

	# Find edges
	fusion_annotation <- gene_file[gene_file[,1] == fusion, , drop = FALSE]
	lower <- 0
	upper <- max(unlist(fusion_annotation["end"]))

	# return list of edges
	multiple_return <- list(lower = lower, upper = upper)
	return(multiple_return)
}

#---------------------------------------------------
#
#   Prepare Annotation Tracks
#
#   Input: Data frame of annotation file contents
#   Output: Return group for track annotation
#
#---------------------------------------------------

prepAnnotation <- function(file){

	raw_names <- as.character(unlist(file[,4]))
	names <- gsub("Exon", "", raw_names)
	count <- nrow(file)
	individuals <- rep(1, count)

	multiple_return <- list(group = rep(names, individuals), count = count)
	return(multiple_return)
}

prepAnnotationDomain <- function(file){

	count <- nrow(file)
	names <- sprintf(paste("Domain", "%d", sep=""), 1:count)

	individuals <- rep(1, count)

	multiple_return <- list(group = rep(names, individuals), count = count)
	return(multiple_return)
}

prepAnnotationTranscripts <- function(file){

	count <- nrow(file)
	names <- sprintf("%d", 1:count)

	individuals <- rep(1, count)

	multiple_return <- list(group = rep(names, individuals), count = count)
	return(multiple_return)
}

#---------------------------------------------------
#
#   Colour Annotation Tracks
#
#   Input: Data Frame
#   Output: Return group for track annotation
#
#---------------------------------------------------

colourTracks <- function(name, count, random, ...){

	ids <- sprintf(paste(name, "%d", sep=""), 1:count)
	colours <- list(...)

	if(random){
		colour_frame <- data.frame(name = ids, colour = randomColour(count))
	} else {
		colour_frame <- data.frame(name = ids, colour = rep(c(colours[[1]], colours[[2]]), count))
	}

	coloured <- paste(colour_frame$name,'="',colour_frame$colour,'"', collapse = ",", sep="")
	multiple_return <- list(colour = coloured, ids = ids)

	return(multiple_return)
}

#---------------------------------------------------
#
#   Allocate Genes based on 1 / 0 from annotation
#
#   Input: Data Frame
#   Output: Return group for track annotation
#
#---------------------------------------------------

geneBoundaries <- function(annotations){

	# Find Gene Boundaries
	gene_allocation <- annotations$genes[,5]
	gene_allocation <- as.character(gene_allocation)
	gene_allocation[gene_allocation == "0"] <- "gene_one"
	gene_allocation[gene_allocation == "1"] <- "gene_two"

	return(gene_allocation)
}

#---------------------------------------------------
#
#   Prepare FST plot
#
#   Input: N/A
#   Output: Start and End of alignment/annotation tracks
#
#---------------------------------------------------

prepare <- function(){

	# Get sanitised user input
	user_input <- userInput()
	results_location <- user_input$results_location
	fusion <- user_input$fusion
	fusion_location <- user_input$fusion_location
	alignment_folder <- user_input$alignment_folder

	# Tell the user we're plotting now
	print(paste("Plotting:", fusion, sep=" "))
	print("------------------------------------------------------")

	# Convert fusino name to Windows friendly format
	fusion_friendly <- gsub(":", "_", fusion)

	# Rename normalised coverage... makes life a bit easier
	old_normalised_alignment = paste(alignment_folder, "Signal.UniqueMultiple.str1.out.bg", sep="/")
	normalised_alignment = paste(alignment_folder, "coverage_rpm.bedgraph", sep="/")

	if(file.exists(old_normalised_alignment)){
		file.rename(old_normalised_alignment, normalised_alignment)
	}

	# Locate alignment and annotation files
	splice_junction_reads <- paste(fusion_location, "splice_junction_reads.bam", sep="/")
	fusion_breakpoint_reads <- paste(fusion_location, "fusion_breakpoint_reads.bam", sep="/")
	gene_location <- paste(results_location, "annotation/gene_boundaries.bed", sep="/")
	transcripts_location <- paste(results_location, "annotation/transcripts.gtf", sep="/")
	protein_location <- paste(results_location, "annotation/protein_domains.bed", sep="/")
	sj_tab_location <- paste(fusion_location, "junctions.txt", sep="/")

	# Put all of these locations in a simple list
	locations <- list(  coverage = normalised_alignment,
						splice_junctions = splice_junction_reads,
						fusion_breakpoints = fusion_breakpoint_reads,
						genes = gene_location,
						transcripts = transcripts_location,
						junctions = sj_tab_location,
						proteins = protein_location)

	# Update the user on progress
	print("Libraries and ancillary files loaded. Creating Tracks.")

	# Decision time, single gene or fusion (does it have a :)?
	if(grepl(":", fusion, fixed=TRUE)){
		is_fusion <- TRUE
	} else {
		is_fusion <- FALSE
	}

	# Load in gene annotation
	gene_file <- as.data.frame(read.table(locations$genes))

	if(is_fusion){
		gene_boundaries <- gene_file[gene_file[,1] == fusion, , drop=FALSE] 
	} else {
		gene_filtered <- unique(gene_file[gene_file$V4 == fusion, , drop=FALSE])
		gene_filtered <- gene_filtered[1,,drop=FALSE]
		gene_boundaries <- gene_filtered[1,]
		gene_boundaries$V1 <- fusion
	}

	colnames(gene_boundaries) <- c("chromosomes", "start", "end", "name", "gene_no", "strand", "start2", "end2", "rgb")

	# Load in the transcripts
	transcript_file <- as.data.frame(read.table(locations$transcripts))

	if(is_fusion){
		transcripts_filtered <- transcript_file[transcript_file[,1] == fusion, , drop=FALSE]
		transcript_boundaries <- transcripts_filtered[,c(1,4,5,7,13,16)]
	} else {
		transcripts_filtered <- transcript_file[grepl(fusion, transcript_file$V1), ]
		transcripts_fusion <- transcripts_filtered[1,1, drop=FALSE]
		transcripts_filtered <- transcripts_filtered[transcripts_filtered$V1 == transcripts_fusion[1,1] & transcripts_filtered$V5 <= gene_boundaries$end & transcripts_filtered$V4 >= gene_boundaries$start, , drop=FALSE]
		transcript_boundaries <- transcripts_filtered[,c(1,4,5,7,13,16)]
		transcript_boundaries$V1 <- fusion
	}

	colnames(transcript_boundaries) <- c("chromosomes", "start", "end", "strand", "transcript","exon")

	# Load in the splice junctions, differentiate between splice junctions and fusion breakpoints

	if(is_fusion){


		junctions_read <- try(read.table(locations$junctions))
		if(class(junctions_read)=='try-error') {
			print("There are no fusion junctions found in this gene pair.")
			print("Terminating plotting for this fusion.")
			return()
		}

		junctions_file <- as.data.frame(junctions_read)
		fusion_junctions <- junctions_file[junctions_file[,2] <= gene_boundaries$end[1] & junctions_file[,3] > gene_boundaries$end[1], , drop=FALSE]

		unique <- fusion_junctions[4]
		multi <- fusion_junctions[5]
		fusion_junctions$support <- unique + multi
		fusion_junction <- fusion_junctions[fusion_junctions$support > 3, , drop=FALSE]

		if(nrow(fusion_junction) < 1){
			stop("Not enough read support for the fusion (< 3 reads supporting), printing terminating.") 
		}

		fusion_frame_name <- c(fusion, fusion)
		fusion_frame_start <- c(fusion_junction[,2], fusion_junction[,3])
		fusion_frame_end <- c(fusion_junction[,2], fusion_junction[,3])
		fusion_frame <- data.frame("chromosomes" = fusion_frame_name, "start" = fusion_frame_start, "end" = fusion_frame_end)

		if(nrow(fusion_junctions) == 0){
			print("There are no fusion junctions found in this gene pair.")
			print("Terminating plotting for this fusion.")
			return()
		}

	} else {
		fusion_frame <- NULL
	}


	# Load in the protein domains
	protein_map_file <- as.data.frame(read.table(locations$proteins))

	if(is_fusion){
		protein_map_filtered <- protein_map_file[protein_map_file[,1] == fusion, , drop=FALSE]
	} else {
		proteins_filtered <- protein_map_file[grepl(fusion, protein_map_file$V1), ]
		protein_fusion <- proteins_filtered[1,1, drop=FALSE]
		protein_map_filtered <- protein_map_file[protein_map_file$V1 == protein_fusion[1,1] & protein_map_file$V3 <= gene_boundaries$end & protein_map_file$V2 >= gene_boundaries$start, , drop=FALSE]
		protein_map_filtered$V1 <- fusion
	}


	colnames(protein_map_filtered) <- c("chromosomes", "start", "end", "domain")

	# Combined overlapping domains of the same name
	domains_all <- as.data.table(protein_map_filtered)
	domains_all[,group := {ir <-  IRanges(start, end); subjectHits(findOverlaps(ir, reduce(ir)))}, by=domain] 

	## Now I group by group and chrom 
	protein_map_merged <- data.frame(domains_all[, list(chromosomes=list(chromosomes), start=min(start),end=max(end), domain=unique(domain)), by=list(group, domain)])
	protein_map_merged <- protein_map_merged[, c("chromosomes","start","end","domain")]
	protein_map_merged$chromosomes <- fusion

	# FINALLY! Let's send the annotations off in style
	annotations <- list(genes = gene_boundaries,
					proteins = protein_map_merged,
					junctions = fusion_frame,
					transcripts = transcript_boundaries)


	create(locations, annotations, results_location, fusion, fusion_friendly, user_input)

}

#---------------------------------------------------
#
#   Setup Plot
#
#   Input: Results of prepare
#   Output: Coloured tracks read for plotTracks
#
#---------------------------------------------------

create <- function(locations, annotations, results_location, fusion, fusion_friendly, user_input){

	user_track_colours <- unlist(strsplit(user_input$track_colours, ","))
	user_track_order <- unlist(strsplit(user_input$track_order, ","))
	track_colours <- paste("#", user_track_colours, sep="")

	if(grepl(":", fusion, fixed=TRUE)){
		is_fusion <- TRUE
	} else {
		is_fusion <- FALSE
	}

	# Prepare Annotation Tracks

	gene_group <- prepAnnotation(annotations$genes)

	protein_group <- prepAnnotationDomain(annotations$proteins)
	protein_id <- prepAnnotation(annotations$proteins)

	transcript_group <- prepAnnotationTranscripts(annotations$transcripts)
	highlight_start <- unname(unlist(annotations$junctions["start"]))
	highlight_end <- unname(unlist(annotations$junctions["end"]))

	# Track Colours - coverage, gene 1, gene 2, domains, transcripts, transcript_bg, transcript_border, junctions, fusions

	# Gene Axis Track
	axis <- GenomeAxisTrack(fontsize = 12, chromosome = fusion, add53 = TRUE, add35 = TRUE)

	# Normalised coverage track
	coverage <- DataTrack(  
		locations$coverage, 
		fontsize = 10, 
		fontsize.title = 12,
		background.title = track_colours[1], 
		name = "Coverage (RPM)", 
		chromosome = fusion, 
		type = c("histogram"), 
		fill = track_colours[1], 
		col.histogram = track_colours[1]
	)

	# Gene Annotation Track
	gene_track <- AnnotationTrack(
		annotations$genes, 
		fontsize = 10, 
		fontsize.title = 12,
		showFeatureId = TRUE, 
		showId = FALSE, 
		chromosome = fusion, 
		groupAnnotation = "id", 
		col = track_colours[3], 
		background.title = track_colours[2], 
		id = gene_group$group,
		name="Genes", 
		stacking="hide"
	)

	feature(gene_track) <- c(geneBoundaries(annotations)) # Add names to gene track
		
	# Splice junction track
	splice_junction_track  <- AlignmentsTrack(  
		locations$splice_junctions, 
		sashimiScore = 5, 
		sashimiNumbers = TRUE, 
		fontsize = 10, 
		fontsize.title = 12,
		chromosome = fusion, 
		background.title = track_colours[8], 
		isPaired = T, 
		col.sashimi = track_colours[8], 
		type = c("sashimi"), 
		size = 0.001, 
		lwd = 2, 
		name = "Split Reads", 
		col.axis = track_colours[8]
	)


	# Fusion Breakpoint track
	fusion_breakpoint_track  <- AlignmentsTrack(
		locations$fusion_breakpoints, 
		fontsize = 10,
		fontsize.title = 12, 
		sashimiScore = 5, 
		sashimiNumbers = TRUE, 
		chromosome = fusion, 
		background.title = track_colours[9], 
		isPaired = T, 
		col.sashimi =  track_colours[9], 
		type=c("sashimi"), 
		size = 0.001, 
		lwd = 2, 
		name = "Split Reads", 
		col.axis = track_colours[9]
	)


	# Overlay both splice junction and fusion breakpoint track
	sashimi_plot <- OverlayTrack(
		trackList = list(
			fusion_breakpoint_track, 
			splice_junction_track, 
			fusion_breakpoint_track
		), 
		background.title =  track_colours[9], 
		name = "Split Reads", 
		col.axis = track_colours[9]
	)

	# Protein Domain Annotation Track
	domain_track <- AnnotationTrack(
		annotations$proteins, 
		fill= track_colours[4], 
		col = track_colours[4], 
		fontsize = 10, 
		fontsize.title = 12,
		groupAnnotation = "id", 
		showId=TRUE, 
		just.group = "below", 
		fontcolor.item = "#000000", 
		background.title = track_colours[4], 
		title.width = 0.5,  
		chromosome = fusion, 
		id = protein_id$group, 
		group = protein_group$group, 
		name="Domains"
	)

	# Transcript Annotation Track
	transcript_track <- GeneRegionTrack(
		annotations$transcripts, 
		chromosome = fusion, 
		fontsize = 10, 
		fontsize.title=12, 
		background.panel = track_colours[6], 
		showExonId=TRUE, 
		showId=TRUE, 
		fill = track_colours[5], 
		col = track_colours[7], 
		just.group = "below", 
		background.title = track_colours[5], 
		title.width = 0.5, 
		name="Transcripts"
	)

	# Create lines over fusion breakpoints

	highlights <- HighlightTrack(
		list(coverage, gene_track, domain_track, transcript_track),
		start = highlight_start, 
		end = highlight_end, 
		col = "transparent", 
		fill = track_colours[10], 
		inBackground = FALSE, 
		chromosome = fusion
	)

	# Create PDF
	pdf_width <- as.integer(user_input$pdf_width)
	pdf_height <- as.integer(user_input$pdf_height)
	tracks <- list(axis, highlights, sashimi_plot)

	pdfIt(pdf_width, pdf_height, fusion, fusion_friendly, results_location, annotations, tracks, user_input)

}

#---------------------------------------------------
#
#   PDF IT!
#
#   Input: Dimensions and results of create()
#   Output: PDF printed into fusion friendly folder
#
#---------------------------------------------------

pdfIt <- function(width, height, fusion, fusion_friendly, results_location, annotations, tracks, user_input){

	user_track_colours <- unlist(strsplit(user_input$track_colours, ","))
	track_colours <- paste("#", user_track_colours, sep="")

	print("Tracks created, printing PDF.")
	sample_name <- user_input$sample_name
	pdf_name <- paste(fusion_friendly, '.pdf', sep="")

	if(is.na(sample_name)){
		pdf_location <- paste(results_location, "plots", pdf_name, sep="/")
	} else {
		pdf_location <- paste(results_location, "plots", sample_name, pdf_name, sep="/")
	}

	# Get dimensions for final track

	track_dimensions <- trackDimensions(annotations$gene, fusion)

	# Create the PDF
	pdf(pdf_location, width=width, height=height)

	ratio <- as.integer(strsplit(user_input$ratio, ",")[[1]])

	print(ratio)

	# Add content
	plot_it <- 'plotTracks(tracks, add=TRUE, sizes = ratio,
				from=0, to=track_dimensions$upper,
				chromosome=fusion, shape = "box", min.height = 5,
				gene_one=track_colours[2], gene_two=track_colours[3])'

	the_plot <- eval(parse(text = plot_it))

	print("PDF created.")

	# Clean up
	garbage <- dev.off()

}

#-----------------------------------------------------------
#
#   R U N
#
#-----------------------------------------------------------

prepare()
