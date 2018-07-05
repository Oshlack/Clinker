'''=========================================================

    L I B R A R I E S

========================================================='''

import setup
import sys
from random import randint

'''=========================================================

    C L A S S E S

========================================================='''

class Fusions:

    def __init__(self):
        self.list = []

    def addFusion(self, fusion):
        self.list.append(fusion)

    def noFusions(self):
        print len(self.list)

class Fusion:

    def __init__(self, gene_name_1, gene_name_2, gene_sequence_1, gene_sequence_2):
        self.gene_name_1 = gene_name_1
        self.gene_name_2 = gene_name_2
        self.gene_sequence_1 = gene_sequence_1
        self.gene_sequence_2 = gene_sequence_2
        self.fused_sequence = gene_sequence_1+gene_sequence_2

class Chromosomes:

    def __init__(self):
        self.index = {}

    def addChromosome(self, name, gene_list):
        try:
            current_gene_list = self.index[name]
            self.index[name] = current_gene_list + gene_list
        except KeyError:
            self.index[name] = gene_list

class Gene:

    def __init__(self, name, start, end):
        self.name = name
        self.start = start
        self.end = end


'''=========================================================

    F U N C T I O N S

========================================================='''

'''---------------------------------------------------------
#
#   Generate a response based upon a triggered error code
#
#   Input:    Error code
#   Output:   Error message, potentially an exit
#
---------------------------------------------------------'''

def error(code, *support):

    '''
    CODE FORMAT:
    First number = where (1 ~ pos, 2 ~ gene, 3 ~ superTrascript)
    Second number = type (1 ~ input, 2 ~ mapping, 3 ~ message)
    '''

    if code == 11:
        print "\nERROR: Column assignment failure. Are you sure you selected the correct columns, are they in the right format?"
        print "Please ensure following format: -pos c4hr1_col breakpoint_1_col chr2_col breakpoint_2_col \n"
        print "==============================================================\n"
        sys.exit(1)

    if code == 22:
        min_gene = support[0]
        min_distance = support[1]
        total = support[2]
        stream = support[3]
        new_line = ''.rjust(9, ' ')

        print "WARNING: a gene (line "+ str(total) +" of fusion input) does not exist in annotation/hg19_ucscGenes.txt based upon breakpoint."
        print new_line+"Closest mapped gene name is '"+min_gene+"' ("+str(min_distance)+" bp "+stream+")"

    if code == 23:
        print "--------------------------------------------------------------\n"
        print "Note: Some superTranscripts were not generated. This could be because of:"
        print "\tA: The breakpoint was not within a gene (this program only deals with these)."
        print "\tB: The superTranscript reference file did not contain an entry for that gene symbol."
        print "\tC: You have identified the wrong columns, or they contain the wrong information, with the -pos argument.\n"
        print "==============================================================\n"

    if code == 32:
        print "WARNING: " + str(support[0]) + " Not found in superTranscriptome"

'''---------------------------------------------------------
#
#   Create a gene list for hg19 or hg38
#
#   Input:    Gene List File
#   Output:   Chromosome object
#
---------------------------------------------------------'''

def geneList(gene_locations):

    first = True
    c_chr = False
    chromosomes = Chromosomes()

    for line in gene_locations:
        gene_line = line.split("\t")
        gene = Gene(gene_line[1], int(gene_line[2]), int(gene_line[3]))

        if gene_line[0].lower() != c_chr:

            if first:
                first = False
            else:
                chromosomes.addChromosome(c_chr, gene_list)

            gene_list = []
            c_chr = gene_line[0].lower()
            gene_list.append(gene)

        else:
            gene_list.append(gene)

    return chromosomes

'''---------------------------------------------------------
#
#   Get chromosome names and breakpoint coord
#
#   Input:      Position argument
#   Output:     Chromosomes and breakpoint coordinates from a
#               fusion finder result
#
---------------------------------------------------------'''

def mapGene(chromosomes, chromosome, breakpoint, total):

    gene_default = "Not Found"
    min_distance = 0
    min_gene = False
    gene_length_min = 9000000000
    found_gene = False

    try:

        for gene in chromosomes.index[chromosome]:

            # For nested genes, choose the smallest gene in which the breakpoint resides
            gene_length = abs(gene.end - gene.start)

            if breakpoint >= gene.start and breakpoint <= gene.end:
                if gene_length < gene_length_min:
                    confirmed_gene = gene.name
                    gene_length_min = gene_length
                    found_gene = True
            else:
                start_distance = abs(gene.start - breakpoint)
                end_distance = abs(gene.end - breakpoint)

                if start_distance < end_distance:
                    if min_gene == False or start_distance < min_distance:
                        min_distance = start_distance
                        min_gene = gene.name
                        where = "upstream"
                elif start_distance > end_distance:
                    if min_gene == False or end_distance < min_distance:
                        min_distance = end_distance
                        min_gene = gene.name
                        where = "downstream"

        if found_gene:
            return confirmed_gene
        else:
            error(22, min_gene, min_distance, total, where)
            return gene_default

    except KeyError:
        return gene_default

'''---------------------------------------------------------
#
#   Get Sequence
#
#   Input:      List of superTranscripts, gene of interest
#   Output:     Sequence
#
---------------------------------------------------------'''

def getSequence(st_genes, gene):

    try:
        sequence = st_genes[gene][0]
    except KeyError:
        sequence = False
        #error(32, gene)

    return sequence


'''---------------------------------------------------------
#
#   Create fused superTranscript
#
#   Input:      Gene 1 and 2
#   Output:     List of fused superTranscripts
#
---------------------------------------------------------'''

def mapSupertranscript(gene_1, gene_2, duplicates, found, not_found, total, st_genes):

    if gene_1 == "Not Found" or gene_2 == "Not Found":
        not_found +=1
        return False, duplicates, found, not_found
    else:
        found += 1

        # Check if this is a duplicate
        try:
            if duplicates[gene_1+":"+gene_2] == True:
                return False, duplicates, found, not_found
        except KeyError:
            duplicates[gene_1+":"+gene_2] = True

        # Check if both genes are in the superTranscriptome
        sequence_1 = getSequence(st_genes, gene_1)
        sequence_2 = getSequence(st_genes, gene_2)

        # We have the gene names and sequences, create a fusion and add it to the list of fusions
        if sequence_1 != False and sequence_2 != False:
            fusion = Fusion(gene_1, gene_2, sequence_1, sequence_2)
            return fusion, duplicates, found, not_found
        else:
            return False, duplicates, found, not_found


'''---------------------------------------------------------
#
#   Get chromosome names and breakpoint coord
#
#   Input:      Position argument
#   Output:     Chromosomes and breakpoint coordinates from a
#               fusion finder result
#
---------------------------------------------------------'''

def breakpointLocation(pos, columns):

    if len(pos) == 2:

        # Map breakpoints to gene symbols
        col_chr_1 = pos[0]
        col_chr_2 = pos[1]

        # Get chromosome number and breakpoint coordinates
        if ";" in columns[int(col_chr_1) - 1]:
            chr_entry_1 = columns[int(col_chr_1) - 1].split(";")
        else:
            chr_entry_1 = columns[int(col_chr_1) - 1].split(":")

        if ";" in columns[int(col_chr_2) - 1]:
            chr_entry_2 = columns[int(col_chr_2) - 1].split(";")
        else:
            chr_entry_2 = columns[int(col_chr_2) - 1].split(":")

        try:
            chr_1 = chr_entry_1[0].strip('"').strip("'").lower()
            chr_2 = chr_entry_2[0].strip('"').strip("'").lower()
            bp_1 = int(chr_entry_1[1])
            bp_2 = int(chr_entry_2[1])
        except:
            error(11)

    else:

        # Map breakpoints to gene symbols
        col_chr_1 = pos[0]
        col_bp_1 = pos[1]
        col_chr_2 = pos[2]
        col_bp_2 = pos[3]

        # Get chromosome number and breakpoint coordinates
        try:
            chr_1 = columns[int(col_chr_1) - 1].strip('"').strip("'").lower()
            bp_1 = int(columns[int(col_bp_1) - 1])
            chr_2 = columns[int(col_chr_2) - 1].strip('"').strip("'").lower()
            bp_2 = int(columns[int(col_bp_2) - 1])
        except:
            error(11)

    # Ensure a common chr format, i.e. chr17
    if "chr" not in chr_1:
        chr_1 = "chr"+chr_1
        chr_2 = "chr"+chr_2

    return chr_1, bp_1, chr_2, bp_2

'''---------------------------------------------------------
#
#   Create exon_boundaries.bed, gene_boundaries.bed
#
#   Input:    Opened Caller File, genes list, annotation path
#
---------------------------------------------------------'''

def createAnnotationFiles(fusions, st_genes, annotation_folder):

    exon_boundaries = open(annotation_folder+'/exon_boundaries.bed','w')
    gene_boundaries =  open(annotation_folder+'/gene_boundaries.bed','w')

    for fusion in fusions.list:

        gene_1_exons = st_genes[fusion.gene_name_1][1].split(",")

        exon_no = 0

        for exons in gene_1_exons:

            exon_no += 1

            exon_range = exons.split("-")
            exon_start = int(exon_range[0]) - 1
            exon_end = int(exon_range[1]) - 1

            R = str(randint(0,255))
            G = str(randint(0,255))
            B = str(randint(0,255))

            RGB = R+","+G+","+B

            exons = fusion.gene_name_1 + ":" + fusion.gene_name_2 + "\t" + str(exon_start) + "\t" + str(exon_end) + "\t" + "Exon"+str(exon_no) + "\t" + "0" + "\t" + "." + "\t" +  str(exon_start) + "\t" +  str(exon_end) + "\t" + RGB + "\n"
            exon_boundaries.write(exons)

        gene_line = fusion.gene_name_1 + ":" + fusion.gene_name_2+"\t0\t"+str(exon_end + 1)+"\t"+fusion.gene_name_1+"\t0\t+\t0\t"+str(exon_end + 1)+"\t255,0,0\n"
        gene_boundaries.write(gene_line)

        gene_1_end = int(exon_end)
        gene_2_start = int(exon_end) + 1

        gene_2_exons = st_genes[fusion.gene_name_2][1].split(",")
        exon_no = 0

        for exons in gene_2_exons:

            exon_no += 1

            exon_range = exons.split("-")
            exon_start = int(gene_1_end) + int(exon_range[0]) - 1
            exon_end = int(gene_1_end) +int(exon_range[1]) - 1

            R = str(randint(0,255))
            G = str(randint(0,255))
            B = str(randint(0,255))

            RGB = R+","+G+","+B

            exons = fusion.gene_name_1 + ":" + fusion.gene_name_2 + "\t" + str(exon_start) + "\t" + str(exon_end) + "\t" + "Exon"+str(exon_no) + "\t" + "1" + "\t" + "." + "\t" +  str(exon_start) + "\t" +  str(exon_end) + "\t" + RGB + "\n"
            exon_boundaries.write(exons)


        gene_line = fusion.gene_name_1 + ":" + fusion.gene_name_2+"\t"+str(gene_2_start)+"\t"+str(exon_end)+"\t"+fusion.gene_name_2+"\t1\t+\t"+str(gene_2_start)+"\t"+str(exon_end)+"\t0,186,255\n"
        gene_boundaries.write(gene_line)


    exon_boundaries.close()
    gene_boundaries.close()


'''---------------------------------------------------------
#
#   Create FASTA file through coordinate / gene lookup
#
#   Input:    Opened Caller File
#   Output:   Column id for fusion genes and class
#
---------------------------------------------------------'''

def createFusionList(fusion_results, pos, gene_list_location, st_genes, header, delimiter):

    if type(fusion_results) == str:
        supplied = True
        fusion_results = fusion_results.strip("[").strip("]").split(",")
    else:
        supplied = False

    fusions = Fusions()
    duplicates = {}

    # Lets index that hg19 or hg38 gene descriptor
    gene_locations = open(gene_list_location, "r")
    chromosomes = geneList(gene_locations)

    # Start tallying warnings
    found = 0
    not_found = 0
    total = 0

    # Error duplicates
    st_errors = {}

    # Print stage
    print "\n==============================================================\n"
    print "Create fusion superTranscriptome:\n"

    for line in fusion_results:

        if not supplied:

            if header:
                header = False
                continue

            if delimiter == "t":
                delim = "\t"
            else:
                delim = ","

            columns = line.split(delim)

            # If two columns selected, not 4, then chromosome and location in same cell. Split by ; and :
            # Assumption that genomic coordinates are in the form of chrX:1932801
            chr_1, bp_1, chr_2, bp_2 = breakpointLocation(pos, columns)

            # Map gene names from breakpoint coordinates
            gene_1 = mapGene(chromosomes, chr_1, bp_1, total)
            gene_2 = mapGene(chromosomes, chr_2, bp_2, total)

        else:
            gene_entry = line.split(":")

            gene_1 = gene_entry[0]
            gene_2 = gene_entry[1]

        # Map gene to superTranscriptome
        fusion, duplicates, found, not_found = mapSupertranscript(gene_1, gene_2, duplicates, found, not_found, total, st_genes)

        # Add found fusion superTranscript
        total += 1

        if fusion == False:
            continue
        else:
            fusions.addFusion(fusion)

    # Print Result

    print "\n--------------------------------------------------------------"
    print "Gene Symbols Mapped: "+str(found), "Not Mapped: "+str(not_found), "Total: "+str(total)
    #print "\n--------------------------------------------------------------"
    #print str(len(st_errors))+" genes did not exist within the ST reference file"

    if not_found > 0:
        error(23)
    else:
        print "\n==============================================================\n"

    # cleanup
    if not supplied:
        gene_locations.close()

    return fusions

'''---------------------------------------------------------
#
#   Create FASTA file through comparing Caller results and superTranscriptome
#
#   Input:    Opened Caller File
#   Output:   Column id for fusion genes and class
#
---------------------------------------------------------'''

# Create a fasta file of the found fusions, save into results folder

def createFusionFasta(fusions, reference_folder, st_genes, competitive):

    fusion_st_fasta = open(reference_folder+'/fst_reference.fasta','w')
    first_line = True

    for fusion in fusions.list:
        if not first_line:
            fusion_st_fasta.write("\n")
        else:
            first_line = False

        fusion_st_fasta.write(">" + fusion.gene_name_1 + ":" + fusion.gene_name_2 + "\n")
        fusion_st_fasta.write(fusion.fused_sequence)


    # Add competitive

    first_competitor = True

    if competitive:
        for key in st_genes:

            try:
                st_seq = st_genes[key][0]

                if first_competitor:
                    fusion_st_fasta.write("\n")
                    first_competitor = False

                fusion_st_fasta.write(">" + key + "\n")
                fusion_st_fasta.write(str(st_seq) + "\n")
            except:
                print st_genes[key]

    fusion_st_fasta.close()
