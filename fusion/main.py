'''

F U S I O N   S U P E R T R A N S C R I P T

Author: Breon Schmidt
Date: 21 Nov 2016

Create a superTranscript (st) of high confidence fusions found through fusion Caller

INPUT: .fasta, Caller.csv
OUTPUT: superTranscript, .bed file

'''

'''=========================================================

	L I B R A R I E S

========================================================='''

import sys
import os.path

import setup
import fusiontools
import annotation

'''=========================================================

    R U N

========================================================='''

user_input = sys.argv

# Get input
fusion_input, destination, pos, gene_list_location, genomic_coordinates, genome_build_location, st_location, annotation_location, protein_location, delimiter, competitive, supplied_fusions, header = setup.programInput(user_input)

# Create a hash table of super transcripts
st_file = open(st_location, 'r')
st_genes = setup.superTranscriptIndex(st_file)

# clean up
st_file.close()

# go through Fusion Finders results
if supplied_fusions == False:
	fusion_results = open(fusion_input,'r')
else:
	fusion_results = supplied_fusions

if "pos" in locals():
	fusions = fusiontools.createFusionList(fusion_results, pos, gene_list_location, st_genes, header, delimiter)
else:
	pos = [1, 2, 3, 4]
	fusions = fusiontools.createFusionList(fusion_results, pos, gene_list_location, st_genes, header, delimiter)


# Create the results folder structure

results_folder = destination
print "Creating output directory at: "+results_folder

if not os.path.exists(results_folder):
    os.makedirs(results_folder)


annotation_folder = results_folder+'/annotation'
reference_folder = results_folder+'/reference'
genome_folder = results_folder+'/genome'
alignment_folder = results_folder+'/alignment'
plots_folder = results_folder+'/plots'

if not os.path.exists(annotation_folder):
	os.makedirs(annotation_folder)
if not os.path.exists(reference_folder):
	os.makedirs(reference_folder)
if not os.path.exists(genome_folder):
	os.makedirs(genome_folder)
if not os.path.exists(plots_folder):
	os.makedirs(plots_folder)
if not os.path.exists(alignment_folder):
	os.makedirs(alignment_folder)

print "Creating fused superTranscriptome and annotation files"

# create a new fasta file containing the fusion supertranscripts, save to results
fusiontools.createFusionFasta(fusions, reference_folder, st_genes, competitive)
fusiontools.createAnnotationFiles(fusions, st_genes, genomic_coordinates, annotation_folder)

# Clean up
if supplied_fusions == False:
	fusion_results.close()

# Create bed files
st_bed_file = open(annotation_location,'r')
st_pbed_file = open(protein_location,'r')
annotation.createAnnotationFile(fusions, st_bed_file, st_pbed_file, annotation_folder)

print "\n...Success!\n\nUse the plot_fst bpipe workflow or IGV to visualise your results."
print "\n==============================================================\n"

#Clean up
st_bed_file.close()
