
'''=========================================================

    C L A S S E S

========================================================='''

class Annotation:

    def __init__(self, annotation): # annotation format taken from https://genome.ucsc.edu/FAQ/FAQformat
        self.chrom = annotation[0]
        self.source = annotation[1]
        self.feature = annotation[2]
        self.chromStart = annotation[3]
        self.chromEnd = annotation[4]
        self.score = annotation[5]
        self.strand = annotation[6]
        self.frame = annotation[7]
        self.group = annotation[8]

    def annotate(self, offset, fusion_name, gene, exon):
        gtf_line = fusion_name+"\t"+self.source+"\t"+self.feature+"\t"+str(int(self.chromStart)+offset)+"\t"+str(int(self.chromEnd)+offset)+"\t"+str(gene)+"\t"+self.strand+"\t"+self.frame+"\t"+self.group+" exon_id "+str(exon)+"\n"
        return gtf_line

class Domain:

    def __init__(self, domain): # annotation format taken from https://genome.ucsc.edu/FAQ/FAQformat

        self.chrom = domain[0]
        self.chromStart = domain[1]
        self.chromEnd = domain[2]
        self.name = domain[3]

    def annotationDomain(self, rgb, fusion_name, offset):

        if int(self.chromStart) < int(self.chromEnd):
           start = int(self.chromStart)
           end = int(self.chromEnd)
        else:
           start = int(self.chromEnd)
           end = int(self.chromStart)

        bed_line = fusion_name+"\t"+str(start+int(offset))+"\t"+str(end+int(offset))+"\t"+self.name+"\n"
        return bed_line


'''=========================================================

    F U N C T I O N S

========================================================='''


'''---------------------------------------------------------
#
#   Create BED file through comparing list of fusions and superTranscriptome BED
#
#   Input:    List of fusions
#   Output:   VOID - creates new file containing annotations
#
---------------------------------------------------------'''

def createAnnotationFile(fusions, st_bed_file, st_pbed_file, annotation_folder):

    # Loop through superTranscriptome, create a dictionary of gene names
    annotations = {}

    for line in st_bed_file:
        annotation = line.strip("\n").split("\t")
        gene_name = annotation[0]

        try:
            annotations[gene_name].append(Annotation(annotation))
        except KeyError:
            annotations[gene_name] = [Annotation(annotation)]

    # Loop through superTranscriptome, create a dictionary of gene names
    domains = {}

    for line in st_pbed_file:
        domain = line.strip("\n").split("\t")
        gene_name = domain[0]

        try:
            domains[gene_name].append(Domain(domain))
        except KeyError:
            domains[gene_name] = [Domain(domain)]

    # Loop through fasta sequence, cross reference against superTranscriptome annotation.
    # Write matches to new file.

    fusion_st_bed_file = open(annotation_folder+'/transcripts.gtf','w')
    fusion_st_pbed_file = open(annotation_folder+'/protein_domains.bed','w')

    for fusion in fusions.list:

        # Check for whether annotation is available
        try:
            annotation_list = annotations[fusion.gene_name_1]
        except KeyError:
            continue

        end_st = 0
        current_transcript = ""

        for annotation in annotation_list:

            current_feature = annotation.feature
            transcript = annotation.group

            if transcript != current_transcript:
                current_transcript = transcript
                exon_count = 1
            else:
                exon_count += 1

            if int(annotation.chromEnd) > end_st:
                end_st = int(annotation.chromEnd)

            gtf_line = annotation.annotate(0, fusion.gene_name_1+":"+fusion.gene_name_2, 0, exon_count);
            fusion_st_bed_file.write(gtf_line)



        # Check for whether annotation is available
        try:
            annotation_list = annotations[fusion.gene_name_2]
        except KeyError:
            continue

        current_transcript = ""

        for annotation in annotation_list:

            current_feature = annotation.feature
            transcript = annotation.group

            if transcript != current_transcript:
                current_transcript = transcript
                exon_count = 1
            else:
                exon_count += 1


            gtf_line = annotation.annotate(end_st, fusion.gene_name_1+":"+fusion.gene_name_2, 1, exon_count);
            fusion_st_bed_file.write(gtf_line)


        # Protein Domains

        try:
            domain_list = domains[fusion.gene_name_1]
        except KeyError:
            continue

        for domain in domain_list:
            pbed_line = domain.annotationDomain("0,255,0", fusion.gene_name_1+":"+fusion.gene_name_2, 0);
            fusion_st_pbed_file.write(pbed_line)

        # Check for whether annotation is available
        try:
            domain_list = domains[fusion.gene_name_2]
        except KeyError:
            continue

        for domain in domain_list:
            pbed_line = domain.annotationDomain("0,255,150", fusion.gene_name_1+":"+fusion.gene_name_2, end_st);
            fusion_st_pbed_file.write(pbed_line)
