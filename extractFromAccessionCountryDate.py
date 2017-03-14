import sys
from Bio import SeqIO, Entrez, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "tpaisie@ufl.edu"
id_list = open("cov_final_accessions.txt")

accessions = []
for line in id_list:
    line = line.rstrip("\n")
    accessions.append(line)

gb_output = open("cov_seqs.gb", "w")
for gb_num in accessions:
    handle = Entrez.efetch(db="nucleotide", id=gb_num, rettype="gb", retmode="text")
    gb_seqs = SeqIO.parse(handle, "gb")
    SeqIO.write(gb_seqs, gb_output, "gb")
handle.close()
gb_output.close()
#genbank file is now downloaded

#extract country and collection date from genbank file
gb_file = open("cov_seqs.gb", "r")
#extract country and collection date from genbank file
#take out accession numbers and country in a text file
save_country = open("cov_country.txt", "w")
for gb_record in SeqIO.parse(gb_file, "genbank") :
    # now do something with the record
        save_country.write(gb_record.name+"\t"),   # print genbank accession number with no new line
        for feat in gb_record.features:
                if feat.type == 'source':
                        source = gb_record.features[0]
                        for qualifiers in source.qualifiers:
                            if qualifiers == 'country':
                                country = source.qualifiers['country']
                                save_country.write(country[0]+"\n"),  #prints the country with no new line
save_country.close()
gb_file.close()

gb_file = open("cov_seqs.gb", "r")
#take out accession numbers and collection date in a text file
# parse genbank file
save_date = open("cov_date.txt", "w")
for gb_record in SeqIO.parse(gb_file, "genbank") :
    # now do something with the record
        save_date.write(gb_record.name+"\t"),   # print genbank accession number with no new line
        for feat in gb_record.features:
                if feat.type == 'source':
                        source = gb_record.features[0]
                        for qualifiers in source.qualifiers:
                                if qualifiers == 'collection_date':
                                        date = source.qualifiers['collection_date']
                                        save_date.write(date[0] + "\n"),  # prints the country with no new line
save_date.close()
gb_file.close()