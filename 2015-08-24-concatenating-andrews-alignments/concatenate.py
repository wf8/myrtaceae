#! /usr/bin/python

import csv
from Bio import Entrez
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from os import listdir
from os.path import isfile, join


print('Concatenating...')
taxa = []
final_sequences = []
alignment_files = [ f for f in listdir('alignments/') if isfile(join('alignments/',f)) ]
records_list = [ list(SeqIO.parse('alignments/' + f, 'fasta')) for f in alignment_files ]

def get_unknowns(length):
    x = ''
    for i in range(length):
        x += '?'
    return x

def parse_binomial(description):
    description = description.replace(' ', '_')
    words = description.split("_")
    if words[1] == 'sp':
        return '_'.join(words[:3])
    else:
        return '_'.join(words[:2])

# get list of all taxa
for records in records_list:
    for record in records:
        taxon = parse_binomial(record.description)
        if taxon not in taxa:
            taxa.append(taxon)
            final_sequences.append('')

# now concatenate
for records in records_list:
    for i, taxon in enumerate(taxa):
        found = False
        for record in records:
            if parse_binomial(record.description) == taxon:
                found = True
                final_sequences[i] += str(record.seq)
                break
        if not found:
            final_sequences[i] += get_unknowns(len(str(records[0].seq)))

# make final seq records
final_records = []
for i in range(len(final_sequences)):
    name = taxa[i].replace(' ', '_') 
    final_records.append(SeqRecord(Seq(final_sequences[i], IUPAC.ambiguous_dna), id=name))

print("Making final phylip-relaxed file...")
SeqIO.write(final_records, "myrtaceae_concatenated.phy", "phylip-relaxed")

print("Done.")
