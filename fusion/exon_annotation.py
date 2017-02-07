'''=========================================================

    L I B R A R I E S

========================================================='''

import sys
from random import randint

'''=========================================================

    R U N

========================================================='''

# Read in the bam file

user_input = sys.argv

annotation_file_location_bed = user_input[1]
exon_annotation = user_input[2]
gene_annotation = user_input[3]

transcripts = open(annotation_file_location_bed, 'r')
#exon_file = open(exon_annotation, 'w')
gene_file = open(gene_annotation, 'w')

new_blocks = []
gene_identifier = []
annotation_name = ""

first_line = True
min_coord = 99999999
max_coord = 0

for transcript in transcripts:

	line = transcript.split("\t")

	gene = line[4]
	block_widths = line[10].split(",")
	block_starts = line[11].split(",")
	del block_starts[-1]
	chrom_start = int(line[1])

	if annotation_name == line[0] or first_line:

		if int(line[1]) < int(min_coord):
			min_coord = line[1]

		if int(line[2]) > int(max_coord):
			max_coord = line[2]

		if first_line:
			annotation_name = line[0]
			first_line = False

		for i in range(len(block_starts)):
			new_blocks.append((chrom_start + int(block_starts[i]), gene))

	else:

		# New genes, filter duplicates
		new_blocks.append((int(max_coord), 1))
		#start = [min_coord]
		#new_blocks = start + new_blocks
		new_blocks = set(new_blocks)
		new_blocks = sorted(new_blocks)
		#new_blocks = map(str, new_blocks)

		k = 0
		tripped = False

		for i in range(1,len(new_blocks)):

			if new_blocks[i-1][1] == str(1) and not tripped:
				k = 1
				tripped = True
				genes = annotation_name.split(":")
				gene_line = annotation_name+"\t0\t"+str(exon_end)+"\t"+genes[0]+"\t0\t+\t0\t"+str(exon_end)+"\t255,0,0\n"
				gene_file.write(gene_line)
				gene_end = exon_end
			else:
				k += 1


			R = str(randint(0,255))
			G = str(randint(0,255))
			B = str(randint(0,255))

			RGB = R+","+G+","+B

			exon_start = new_blocks[i-1][0]
			exon_end = new_blocks[i][0]
			gene_identifier = new_blocks[i-1][1]

			exons = annotation_name + "\t" + str(exon_start) + "\t" + str(exon_end) + "\t" + "Exon"+str(k) + "\t" + str(gene_identifier) + "\t" + "." + "\t" +  str(exon_start) + "\t" +  str(exon_end) + "\t" + RGB + "\n"
			#exon_file.write(exons)


		#  finish gene line
		gene_line = annotation_name+"\t"+str(gene_end)+"\t"+str(exon_end)+"\t"+genes[1]+"\t1\t+\t0\t"+str(exon_end)+"\t0,150,255\n"
		gene_file.write(gene_line)

		# reset
		min_coord = 99999999
		max_coord = 0
		new_blocks = []

		for block_start in block_starts:
			new_blocks.append((chrom_start + int(block_start), gene))

		annotation_name = line[0]


transcripts.close()
#exon_file.close()
gene_file.close()
