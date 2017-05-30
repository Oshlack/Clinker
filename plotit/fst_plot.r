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
        hex_code <- sample(c(0:9,letters[1:6]), 6, replace=TRUE)
        colour_storage[i] <- paste(hex_code, collapse="", sep="")
        colour_storage[i] <- paste("#",colour_storage[i], sep="")
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

    #results_location <- "Path to results folder"
    #fusion <- as.character("Your fusion name")

    # Return multiple variables
    multiple_return <- list(results_location = results_location, fusion = fusion, fusion_location = fusion_location, pdf_width = pdf_width, pdf_height = pdf_height, ratio = ratio, alignment_folder = alignment_folder, sample_name = sample_name)
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
    upper <- max(unlist(fusion_annotation["V3"]))

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

geneBoundaries <- function(files){

    # Find Gene Boundaries
    gene_allocation <- files$genes[,5]
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

    print(paste("Plotting:", fusion, sep=" "))
    print("------------------------------------------------------")

    # Create fusion folder in results
    fusion_friendly <- gsub(":", "_", fusion)

    # Locate alignment and annotation files
    alignment_location <- paste(alignment_folder, "/coverage_rpm.bedgraph", sep="/")
    alignment_location_filtered <- paste(fusion_location, "reads.bam", sep="/")
    splice_reads_location <- paste(fusion_location, "split_reads.bam", sep="/")
    gene_location <- paste(results_location, "annotation/gene_boundaries.bed", sep="/")
    exon_location <- paste(results_location, "annotation/exon_boundaries.bed", sep="/")
    transcripts_location <- paste(results_location, "annotation/transcripts.gtf", sep="/")
    protein_location <- paste(results_location, "annotation/protein_domains.bed", sep="/")
    sj_tab_location <- paste(fusion_location, "junctions.txt", sep="/")

    print("Libraries and ancillary files loaded. Creating Tracks.")

    locations <- list(  alignment = alignment_location,
                        alignment_filtered = alignment_location_filtered,
                        splice = splice_reads_location,
                        genes = gene_location,
                        exons = exon_location,
                        transcripts = transcripts_location,
                        junctions = sj_tab_location,
                        proteins = protein_location)

    # Load in alignment and annotation files
    gene_file <- as.data.frame(read.table(locations$genes))
    gene_filtered <- gene_file[gene_file[,1] == fusion, , drop=FALSE]
    breakpoint <- gene_filtered[1,3]

    exon_file <- as.data.frame(read.table(locations$exons))
    exon_filtered <- exon_file[exon_file[,1] == fusion, , drop=FALSE]
    colnames(exon_filtered) <- c("chromosomes", "start", "end", "exon", "gene", "strandx", "thickStart", "thickEnd", "itemRGB")

    transcript_file <- as.data.frame(read.table(locations$transcripts))
    transcripts_filtered <- transcript_file[transcript_file[,1] == fusion, , drop=FALSE]
    transcript_booundaries <- transcript_file[,c(1,4,5)]
    colnames(transcript_booundaries) <- c("chromosomes", "start", "end")

    junctions_read <- try(read.table(locations$junctions))
    if(class(junctions_read)=='try-error') {
        print("There are no fusion junctions found in this gene pair.")
        print("Terminating plotting for this fusion.")
        return()
    }

    junctions_file <- as.data.frame(junctions_read)
    fusion_junctions <- junctions_file[junctions_file[,2] <= breakpoint & junctions_file[,3] > breakpoint, , drop=FALSE]

    if(nrow(fusion_junctions) == 0){
        print("There are no fusion junctions found in this gene pair.")
        print("Terminating plotting for this fusion.")
        return()
    }

    unique <- fusion_junctions[4]
    multi <- fusion_junctions[5]
    fusion_junctions$support <- unique + multi

    #fusion_junction <- fusion_junctions[fusion_junctions$support == max(fusion_junctions$support), , drop=FALSE]
    print(fusion_junctions$support)
    fusion_junction <- fusion_junctions[fusion_junctions$support > 3, , drop=FALSE]

    fusion_frame_name <- c(fusion, fusion)
    fusion_frame_start <- c(fusion_junction[,2], fusion_junction[,3])
    fusion_frame_end <- c(fusion_junction[,2], fusion_junction[,3])
    fusion_frame <- data.frame("chromosomes" = fusion_frame_name, "start" = fusion_frame_start, "end" = fusion_frame_end)

    protein_map_file <- as.data.frame(read.table(locations$proteins))
    protein_map_filtered <- protein_map_file[protein_map_file[,1] == fusion, , drop=FALSE]
    colnames(protein_map_filtered) <- c("chromosomes", "start", "end", "domain")


    files <- list(  genes = gene_file,
                    exons = exon_filtered,
                    proteins = protein_map_filtered,
                    junctions = fusion_frame,
                    transcripts = transcript_booundaries)

    create(locations, files, results_location, fusion, fusion_friendly, user_input)

}

#---------------------------------------------------
#
#   Setup Plot
#
#   Input: Results of prepare
#   Output: Coloured tracks read for plotTracks
#
#---------------------------------------------------

create <- function(locations, files, results_location, fusion, fusion_friendly, user_input){

    # Prepare Annotation Tracks
    gene_group <- prepAnnotation(files$genes)
    exon_group <- prepAnnotation(files$exons)
    protein_group <- prepAnnotationDomain(files$proteins)
    protein_id <- prepAnnotation(files$proteins)

    highlight_start <- unname(unlist(files$junctions["start"]))
    highlight_end <- unname(unlist(files$junctions["end"]))

    # Splice junction tracks
    splice_junction_track  <- AlignmentsTrack(locations$alignment_filtered, sashimiScore = 10, fontsize = 10, chromosome = fusion, background.title = "#6b98d7", isPaired = T, col.sashimi = "#D7D4E4", type=c("sashimi"), size = 0.001, lwd = 2, name = " ")
    split_read_junction_track  <- AlignmentsTrack(locations$splice, fontsize = 10, sashimiScore = 2, chromosome = fusion, background.title = "#6b98d7", isPaired = T, col.sashimi = "#6e65ad", type=c("sashimi"), size = 0.001, lwd = 2, name = " ")
    sashimi_plot <- OverlayTrack(trackList = list(split_read_junction_track, splice_junction_track, split_read_junction_track), background.title = "#6e65ad", name = " ")

    # Normalised coverage track
    coverage <- DataTrack(locations$alignment, fontsize = 10, background.title = "#6e65ad", name = "Coverage", chromosome = fusion, type=c("histogram"), fill="#6e65ad", col.histogram = "#6e65ad")

    # Annotation tracks

    gene_track <- AnnotationTrack(locations$genes, fontsize = 12, showFeatureId = TRUE, chromosome = fusion, col="#2b749a", background.title = "#3983AA", group = gene_group$group, name="Genes", stacking="hide")
    domain_track <- AnnotationTrack(files$proteins, fill="#f05f3b", col = "#f05f3b", fontsize = 12, groupAnnotation = "id", showId=TRUE, just.group = "below", fontcolor.item = "#000000", background.title = "#f05f3b", title.width = 0.1,  chromosome = fusion, id= protein_id$group, group = protein_group$group, name="Domains")
    transcript_track <- GeneRegionTrack(locations$transcripts, fontcolor="#215A4A", fontsize = 10, fontsize.title=12, background.panel = "#f2ffe4", showExonId=TRUE, showId=TRUE, fill = "#a1d16e", col = "#215A4A", just.group = "below", background.title = "#a1d16e", title.width = 0.1, name="Transcripts")
    exon_track <- AnnotationTrack(files$exons, fontsize = 12, fontsize.group = 8, showFeatureId = TRUE, featureAnnotation = "group", size = 0.1, just.group = "left", fontcolor.item = "#000000", col="black", background.title = "#1d3c66", title.width = 0.1, chromosome = fusion, group = exon_group$group, name="Exons", size = 0.001, stacking = "dense", mergeGroups = FALSE)
    gene_axis_track <- GenomeAxisTrack(fontsize = 12, chromosome = fusion, add53 = TRUE, add35 = TRUE)
    feature(gene_track) <- c(geneBoundaries(files))
    highlights <- HighlightTrack(list(coverage, gene_track, domain_track, transcript_track), start = highlight_start, end = highlight_end, col = "transparent", fill = "#ff6d6d", inBackground = FALSE, chromosome = fusion)

    # Create PDF
    pdf_width <- as.integer(user_input$pdf_width)
    pdf_height <- as.integer(user_input$pdf_height)

    tracks <- list(gene_axis_track, highlights, sashimi_plot)

    pdfIt(pdf_width, pdf_height, fusion, fusion_friendly, results_location, files, tracks, user_input)

}

#---------------------------------------------------
#
#   PDF IT!
#
#   Input: Dimensions and results of create()
#   Output: PDF printed into fusion friendly folder
#
#---------------------------------------------------

pdfIt <- function(width, height, fusion, fusion_friendly, results_location, files, tracks, user_input){

    print("Tracks created, printing PDF.")
    sample_name <- user_input$sample_name

    pdf_name <- paste(fusion_friendly, '.pdf', sep="")

    if(is.na(sample_name)){
        pdf_location <- paste(results_location, "plots", pdf_name, sep="/")
    } else {
        pdf_location <- paste(results_location, "plots", sample_name, pdf_name, sep="/")
    }

    # Get dimensions for final track
    track_dimensions <- trackDimensions(files$gene, fusion)

    # Create the PDF
    pdf(pdf_location, width=width, height=height)
    ratio <- as.integer(strsplit(user_input$ratio, ",")[[1]])

    # Add content
    plot_it <- 'plotTracks(tracks, add=TRUE, sizes = ratio,
                from=track_dimensions$lower, to=track_dimensions$upper,
                chromosome=fusion, shape = "box", min.height = 5,
                gene_one="#3983AA", gene_two="#2b749a")'

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
