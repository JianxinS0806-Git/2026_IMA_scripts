#!/usr/bin/env python

import os
import gffutils


def extract_gene_positions(gff_file):

    current_dir = os.getcwd()
    input_file_path = os.path.join(current_dir, gff_file)
    db = gffutils.create_db(input_file_path, dbfn=':memory:', force=True, keep_order=True, merge_strategy="merge")

    output_file = os.path.splitext(gff_file)[0] + "_gene_info.txt"
    with open(output_file, 'w') as f_out:
        title = "Gene ID" + "\t" + "Contig" + "\t" + "Start" + "\t" + "End" + "\n"
        f_out.write(title)
        for gene in db.features_of_type('gene'):
            gene_id = gene.id
            contig = gene.seqid
            start = gene.start
            end = gene.end
            gene_info = gene_id + "\t" + contig + "\t" + str(start) + "\t" + str(end) + "\n"
            f_out.write(gene_info)

    return output_file

# Get GFF file name
gff_file = input("Please enter the GFF file name: ")

output_file = extract_gene_positions(gff_file)
print("Gene location information has been successfully extracted and written to {}.".format(output_file))
