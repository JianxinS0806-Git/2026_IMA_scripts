#!/usr/bin/env python

from Bio import SeqIO
import numpy as np
import argparse

def calculate_gc_content(cds):
    # Initialize counters
    gc1, gc2, gc3 = 0, 0, 0
    a3, t3, g3, c3 = 0, 0, 0, 0
    
    length = len(cds)
    count = length // 3
    
    for i in range(count):
        codon = cds[i*3:(i+1)*3]
        if len(codon) == 3:
            if codon[0] in 'GC':
                gc1 += 1
            if codon[1] in 'GC':
                gc2 += 1
            if codon[2] in 'GC':
                gc3 += 1
            if codon[2] == 'A':
                a3 += 1
            if codon[2] == 'T':
                t3 += 1
            if codon[2] == 'G':
                g3 += 1
            if codon[2] == 'C':
                c3 += 1
    
    # Calculate ratios
    gc1_ratio = gc1 / count if count > 0 else 0
    gc2_ratio = gc2 / count if count > 0 else 0
    gc3_ratio = gc3 / count if count > 0 else 0
    a3_ratio = a3 / count if count > 0 else 0
    t3_ratio = t3 / count if count > 0 else 0
    g3_ratio = g3 / count if count > 0 else 0
    c3_ratio = c3 / count if count > 0 else 0
    
    gc12_ratio = (gc1_ratio + gc2_ratio) / 2
    at3_ratio = a3_ratio / (a3_ratio + t3_ratio) if (a3_ratio + t3_ratio) > 0 else 0
    gc3_ratio_normalized = g3_ratio / (g3_ratio + c3_ratio) if (g3_ratio + c3_ratio) > 0 else 0
    
    return gc1_ratio, gc2_ratio, gc3_ratio, a3_ratio, t3_ratio, g3_ratio, c3_ratio, gc12_ratio, at3_ratio, gc3_ratio_normalized

def main(input_file, output_file):
    # Read sequences from the input file
    with open(output_file, 'w') as out_f:
        out_f.write("Gene\tGC1\tGC2\tGC3\tA3\tT3\tG3\tC3\tGC12\tA3/(A3+T3)\tG3/(G3+C3)\n")
        
        for record in SeqIO.parse(input_file, "fasta"):
            gene_name = record.id
            cds = str(record.seq)
            
            # Ignore the trailing nucleotides if the length is not a multiple of 3
            if len(cds) % 3 != 0:
                cds = cds[:-(len(cds) % 3)]
            
            results = calculate_gc_content(cds)
            
            out_f.write(f"{gene_name}\t" + "\t".join(map(str, results)) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate GC content for CDS sequences.')
    parser.add_argument("-i", "--input", required=True, help='Input FASTA file containing CDS sequences.')
    parser.add_argument("-o", "--output", required=True, help="Output TSV file to save the results.")

    args = parser.parse_args()
    main(args.input, args.output)
