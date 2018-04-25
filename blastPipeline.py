#running blast over the internet
#query seq already in fasta format, open file and read in this record as a string
import sys
from Bio import SeqIO, Entrez, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML

fasta_string = open(sys.argv[1]).read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string, hitlist_size=100, entrez_query="9500:12000[slen]")
#blast result gives top 100 hits


#saving local copy of the xml output file
save_file = open("my_blast.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()

result_handle = open("my_blast.xml")
#have the handle, read to parse the .xml file
#have multiple Blast results, use parse not read
blast_qresult = SearchIO.read("my_blast.xml", "blast-xml")
#Hit objects represent all query results from a single database entry
#They are the second-level container in the Bio.SearchIO object hierarchy
#Youve seen that they are contained by QueryResult objects
#but they themselves contain HSP objects
save_text = open("blast_accession.txt", "w")
for hit in blast_qresult:
    save_text.write(hit.accession+"\n")
save_text.close()

#from blast_accession.txt use eFetch to download a fasta file of all the accessions
#first uploading the list of IDs using EPost
#can refer to the long list of IDs and download the associated data with EFetch
Entrez.email = input(str(">"))  #enter email here
id_list = open("blast_accession.txt")

accessions = []
for line in id_list:
    line = line.strip()
    accessions.append(line)

with open(sys.argv[2]+".gb", "w") as gb_output:
    for gb_num in accessions:
        with Entrez.efetch(db="nucleotide", id=gb_num, rettype="gb", retmode="text") as handle:
            gb_output.write(handle.read())
            handle.close()

with open(sys.argv[2]+".fasta", "w") as fasta_output:
    for num in accessions:
        with Entrez.efetch(db="nucleotide", id=num, rettype="fasta") as handle:
            fasta_output.write(handle.read())
            handle.close()

#genbank and fasta file are now downloaded

#extract country and collection date from genbank file
gb_file = open(sys.argv[2]+".gb", "r")
#extract country and collection date from genbank file
#take out accession numbers and country in a text file
save_country = open(sys.argv[2]+"_country.txt", "w")
for gb_record in SeqIO.parse(gb_file, "genbank") :
    # now do something with the record
        save_country.write(gb_record.name+" "),   # print genbank accession number with no new line
        for feat in gb_record.features:
                if feat.type == 'source':
                        source = gb_record.features[0]
                        for qualifiers in source.qualifiers:
                            if qualifiers == 'country':
                                country = source.qualifiers['country']
                                save_country.write(country[0]+"\n"),  #prints the country with no new line
save_country.close()
gb_file.close()

gb_file = open(sys.argv[2]+".gb", "r")
#take out accession numbers and collection date in a text file
# parse genbank file
save_date = open(sys.argv[2]+"_date.txt", "w")
for gb_record in SeqIO.parse(gb_file, "genbank") :
    # now do something with the record
        save_date.write(gb_record.name+" "),   # print genbank accession number with no new line
        for feat in gb_record.features:
                if feat.type == 'source':
                        source = gb_record.features[0]
                        for qualifiers in source.qualifiers:
                                if qualifiers == 'collection_date':
                                        date = source.qualifiers['collection_date']
                                        save_date.write(date[0] + "\n"),  # prints the country with no new line
save_date.close()
gb_file.close()

#change the fasta file sequence names to just Genbank accession numbers
#grabbing the file and the names
fasta_file = open(sys.argv[2]+".fasta", "r")
new_fasta = open(sys.argv[2]+"_new_names.fasta", "w")
new_ids = open("blast_accession.txt", "r")

for f in fasta_file.readlines():
    if f.__contains__(">"):
        new_fasta.write(">" + new_ids.readline())
    else:
        new_fasta.write(f)

new_fasta.close()
fasta_file.close()
new_ids.close()

# align fasta file using muscle
from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO

#specify the location of your muscle file
muscle_exe = r"/Users/taylorpaisie/bin/muscle"

input_seqs = sys.argv[2]+"_new_names.fasta"
output_aln = sys.argv[2]+"_aln.fasta"

# tell MUSCLE to read in this FASTA file, and write the alignment to an output file
def align_v1 (Fasta):
    muscle_cline = MuscleCommandline(muscle_exe, input=Fasta, out=output_aln)
    stdout, stderr = muscle_cline()
    MultipleSeqAlignment = AlignIO.read(output_aln, "fasta")
    print(MultipleSeqAlignment)

align_v1(input_seqs)

# aligned fasta file will be made


#to run: python seqBlastPipeline.py FASTAFILE PREFIXFORFILES

