#!/usr/bin/env python

import os

# Extract sequence IDs from the processed protein FASTA file

def extract_sequence_ids(fasta_file):
    sequence_ids = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                sequence_id = line.strip()[1:-3]
                sequence_ids.add(sequence_id)
    return sequence_ids

# Filtering extracted gene information from the GFF file

def filter_gene_info(gene_info_file, sequence_ids):
    filtered_genes = []
    with open(gene_info_file, 'r') as f:
        header = next(f)
        filtered_genes.append(header.strip())
        for line in f:
            gene_id = line.split('\t')[0]
            if gene_id in sequence_ids:
                filtered_genes.append(line.strip())
    return filtered_genes

# Save new file

def save_filtered_gene_info(filtered_genes, output_file):
    with open(output_file, 'w') as f:
        for gene in filtered_genes:
            f.write(gene + '\n')

fasta_file_input = input("Please enter the path or filename of the predicted protein FASTA file: ")
gene_info_file_input = input("Please enter the path or filename of the gene information TXT file extracted from the GFF file: ")

if os.path.isabs(fasta_file_input):
    fasta_file = fasta_file_input
else:
    fasta_file = os.path.abspath(fasta_file_input)

if os.path.isabs(gene_info_file_input):
    gene_info_file = gene_info_file_input
else:
    gene_info_file = os.path.abspath(gene_info_file_input)

sequence_ids = extract_sequence_ids(fasta_file)

filtered_genes = filter_gene_info(gene_info_file, sequence_ids)

output_filename = os.path.splitext(os.path.basename(gene_info_file))[0] + "_filter.txt"
output_file = os.path.join(os.getcwd(), output_filename)

save_filtered_gene_info(filtered_genes, output_file)

print("The filtered gene information has been successfully saved to {}.".format(output_file))
