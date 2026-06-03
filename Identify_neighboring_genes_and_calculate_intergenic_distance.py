#!/usr/bin/env python

import argparse

def main(input_file, output_file):
    gene_positions = {}
    with open(input_file, 'r') as f:
        next(f)  
        for line in f:
            gene_id, contig, start, end = line.strip().split("\t")
            gene_positions[gene_id] = (contig, int(start), int(end))

    gene_distances = {}
    for gene_id, (contig, start, end) in gene_positions.items():
        upstream_gene = None
        upstream_distance = None
        downstream_gene = None
        downstream_distance = None
        
        for other_gene_id, (other_contig, other_start, other_end) in gene_positions.items():
            if contig == other_contig and gene_id != other_gene_id:
                if other_end < start:
                    distance = start - other_end
                    if upstream_distance is None or distance < upstream_distance:
                        upstream_gene = other_gene_id
                        upstream_distance = distance
                elif other_start > end:
                    distance = other_start - end
                    if downstream_distance is None or distance < downstream_distance:
                        downstream_gene = other_gene_id
                        downstream_distance = distance
        
        if upstream_gene is None:
            upstream_gene_data = (None, None, None)
        else:
            upstream_gene_data = (upstream_gene, *gene_positions[upstream_gene][1:])
        
        if downstream_gene is None:
            downstream_gene_data = (None, None, None)
        else:
            downstream_gene_data = (downstream_gene, *gene_positions[downstream_gene][1:])
        
        gene_distances[gene_id] = (start, end, upstream_gene_data, upstream_distance, downstream_gene_data, downstream_distance)

    with open(output_file, "w") as f:
        f.write("Gene ID\tStart\tEnd\tUpstream Gene ID\tUpstream Start\tUpstream End\tUpstream Distance\tDownstream Gene ID\tDownstream Start\tDownstream End\tDownstream Distance\n")
        for gene_id, (start, end, upstream_gene_data, upstream_distance, downstream_gene_data, downstream_distance) in gene_distances.items():
            upstream_gene_id, upstream_start, upstream_end = upstream_gene_data
            downstream_gene_id, downstream_start, downstream_end = downstream_gene_data
            upstream_distance = upstream_distance-1 if upstream_distance is not None else "None"
            downstream_distance = downstream_distance-1 if downstream_distance is not None else "None"
            f.write(f"{gene_id}\t{start}\t{end}\t{upstream_gene_id}\t{upstream_start}\t{upstream_end}\t{upstream_distance}\t{downstream_gene_id}\t{downstream_start}\t{downstream_end}\t{downstream_distance}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate distances between genes and their upstream and downstream neighbors")
    parser.add_argument("input_file", help="Please enter path of the filtered gene information TXT file.")
    parser.add_argument("output_file", help="The output file contains each gene's upstream and downstream neighboring genes and their respective distances.")
    args = parser.parse_args()
    
    main(args.input_file, args.output_file)
