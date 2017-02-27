'''=========================================================

	L I B R A R I E S

========================================================='''

import general

'''=========================================================

	F U N C T I O N S

========================================================='''

def programInput(user_input):

	print "\n\n====================================================\n\n"
	print "\tFusion Super Transcript Generator\n"
	print "\tA fusion visualiser."
	print "\n\n====================================================\n\n"

	for argument in range(1, len(user_input)):

		# Command identifier?
		if "-" == user_input[argument][0]:

			if "-in" == user_input[argument]:
				command = 1
			elif "-out" == user_input[argument]:
				command = 2
			elif "-pos" == user_input[argument]:
				pos = []
				command = 3
			elif "-st" == user_input[argument]:
				command = 4
			elif "-del" == user_input[argument]:
				command = 5
			elif "-genome" == user_input[argument]:
				command = 6
			elif "-competitive" == user_input[argument]:
				command = 7
			elif "-fusionlist" == user_input[argument]:
				command = 8
			else:
				print "ERROR: " + user_input[argument] + " is not a valid argument. Please use the following options:\n"
				print "-in (path/to/fusions OR custom, space delimited list)"
				print "-out (path/to/results_folder)"
				print "-pos (1 2 3 4 - columns for chr1 start chr2 start for your input -in)"
				print "-st (path/to/st/resources, default 'resources')\n"
				sys.exit(1)

		else:

			if command == 1:
				fusion_input = user_input[argument]
			elif command == 2:
				destination = user_input[argument]
			elif command == 3:
				pos.append(user_input[argument])
			elif command == 4:
				resources = user_input[argument]
			elif command == 5:
				delimiter = user_input[argument]
			elif command == 6:
				genome_build = user_input[argument]
			elif command == 7:
				if user_input[argument] == "true":
					competitive = True
				elif user_input[argument] == "false":
					competitive = False
				else:
					print "ERROR: Please ensure that -competitive is set to either true or false"
					sys.exit(1)
			elif command == 8:
				if len(user_input[argument]) > 0:
					supplied_fusions = user_input[argument]
				else:
					supplied_fusions = False

	# If resources does not exist, assume native structure
	if "resources" not in locals():
		resources = "../resources"
		st_location = resources+'/hg38_genCode24_st.fasta'
		annotation_location = resources+'/hg38_genCode24_st-sorted-exons.gtf'
		protein_location = resources+'/hg38_genCode24_pfam-map.bed'

	if "genome_build" in locals():
		if genome_build == "19":
			genome_build_location = "hg19_ucscGenes.txt"
		elif genome_build == "38":
			genome_build_location = "hg38_genCode24.txt"
		else:
			print "ERROR: Genome build can be either 38 or 19, please select appropriately."
			sys.exit(1)

		gene_list_location = resources+'/'+genome_build_location

	if "delimiter" not in locals():
		delimiter = "c"

	if "supplied_fusions" not in locals():
		supplied_fusions = False

	return fusion_input, destination, pos, gene_list_location, genome_build_location, st_location, annotation_location, protein_location, delimiter, competitive, supplied_fusions

'''---------------------------------------------------------
#
#   Create a hash table of a superTranscriptome
#
#   Input:	An opened, readable file
#   Output:   A hash table with gene names as the key,
#			 values are the sequence
#
---------------------------------------------------------'''

def superTranscriptIndex(st_file):

	st_genes = {}
	first_line = True

	# don't use readlines, commits whole file to memory
	for line in st_file:

		# if the first characters is >, gene name, else, sequence.
		first_char = line[0]
		if first_char == ">":

			# Assign the previous gene with a concenated sequence (if line breaks present)
			if not first_line:
				full_sequence = "".join(sequence)
				st_genes[gene_name] = (full_sequence, exon_coordinates)
			else:
				first_line = False

			# remove the gene indicator, ">", spaces and breaks
			line_parts = line.split(" ")
			gene_name = general.stripStringChars(line_parts[0], ["\n", ">"])
			exon_coordinates = line_parts[1].strip("segs:")
			st_genes[gene_name] = ""

			# Start collating a new sequence
			sequence = []

		else:
			# remove any breaks or white space from sequence
			sequence_line = general.stripStringChars(line, ["\n", " "])
			sequence.append(sequence_line)

	return st_genes
