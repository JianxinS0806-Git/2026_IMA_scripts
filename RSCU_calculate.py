#!/usr/bin/env python

import argparse
from Bio import SeqIO
from collections import Counter, defaultdict

# Calculate RSCU

def calculate_rscu(codon_counts, codon_table):
    rscu_values = {}
    for amino_acid, codon_list in codon_table.items():
        total_count = sum(codon_counts[codon] for codon in codon_list)
        num_codons = len(codon_list)
        for codon in codon_list:
            observed_count = codon_counts[codon]
            expected_count = total_count / num_codons if num_codons > 0 else 0
            rscu = observed_count / expected_count if expected_count > 0 else 0
            rscu_values[codon] = rscu
    return rscu_values

# Process sequence

def process_sequence(sequence):
    # Remove start codon (first three nucleotides) and any stop codons (last three nucleotides if they form a stop codon)
    sequence = sequence[3:]  # remove start codon
    if len(sequence) % 3 != 0:
        sequence = sequence[:-(len(sequence) % 3)]  # remove incomplete codon at the end
    stop_codons = {'TAA', 'TAG', 'TGA'}
    if sequence[-3:] in stop_codons:
        sequence = sequence[:-3]
    return sequence

def main(input_file, output_file):
    codon_table = {
        'F': ['TTT', 'TTC'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'I': ['ATT', 'ATC', 'ATA'],
        'M': ['ATG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'Y': ['TAT', 'TAC'],
        'H': ['CAT', 'CAC'],
        'Q': ['CAA', 'CAG'],
        'N': ['AAT', 'AAC'],
        'K': ['AAA', 'AAG'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'C': ['TGT', 'TGC'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'W': ['TGG']
    }
    all_codons = [codon for codons in codon_table.values() for codon in codons]

    species_codon_counts = Counter()
    gene_codon_counts = defaultdict(Counter)
    
    for record in SeqIO.parse(input_file, "fasta"):
        gene_name = record.id
        sequence = str(record.seq).upper()
        processed_sequence = process_sequence(sequence)
        
        codons = [processed_sequence[i:i+3] for i in range(0, len(processed_sequence), 3)]
        gene_codon_counts[gene_name].update(codons)
        species_codon_counts.update(codons)
    
    species_rscu = calculate_rscu(species_codon_counts, codon_table)
    gene_rscu = {gene: calculate_rscu(codon_counts, codon_table) for gene, codon_counts in gene_codon_counts.items()}
    
    with open(output_file, 'w') as f:
        header = ["gene"] + all_codons
        f.write("\t".join(header) + "\n")
        
        species_line = ["species"] + [str(species_rscu.get(codon, 0)) for codon in all_codons]
        f.write("\t".join(species_line) + "\n")
        
        for gene, rscu_values in gene_rscu.items():
            gene_line = [gene] + [str(rscu_values.get(codon, 0)) for codon in all_codons]
            f.write("\t".join(gene_line) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RSCU values from a FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output file for RSCU values")
    args = parser.parse_args()
    
    main(args.input, args.output)
