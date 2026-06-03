#!/usr/bin/env python

import os

# Read sequences

def read_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        seq = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    sequences.append(''.join(seq))
                    seq = []
            else:
                seq.append(line)
        if seq:
            sequences.append(''.join(seq))
    return sequences

# Create sequence matrix

def create_sequence_matrix(directory):
    matrix = []
    files = sorted([os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.fasta')])

    for file in files:
        sequences = read_fasta(file)
        matrix.append(sequences)

    return matrix

# Concatenate sequences column-wise in order

def concatenate_sequences(sequence_matrix):
    num_genes = len(sequence_matrix[0])
    concatenated_sequences = ['' for _ in range(num_genes)]
    
    for row in sequence_matrix:
        for i in range(num_genes):
            concatenated_sequences[i] += row[i]
    
    return concatenated_sequences

def write_concatenated_sequences_to_fasta(concatenated_sequences, output_file):
    with open(output_file, 'w') as f:
        for i, sequence in enumerate(concatenated_sequences):
            f.write(f">ID{i+1}\n")
            f.write(f"{sequence}\n")

# Directory containing MSA files for all single-copy orthologous genes, each file corresponds to one aliged single-copy ortholog.
directory_path = "/home/jxs/Working/C_result/protein/trimal/"

sequence_matrix = create_sequence_matrix(directory_path)

concatenated_sequences = concatenate_sequences(sequence_matrix)

output_file = "concatenated_sequences.fasta"
write_concatenated_sequences_to_fasta(concatenated_sequences, output_file)

print(f"Concatenated sequences have been written to {output_file}")

