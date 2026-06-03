This repository contains the scripts used in the analyses presented in the study "Genus-wide comparative genomics of Colletotrichum reveals evolutionary conservation and divergent ecological adaptations".

1. Calculate intergenic distances

Extract_gene_information_from_GFF.py
Extracts all predicted gene information from an annotated GFF file, including gene IDs, contig names, and the start and end coordinates of each gene.

Gene_information_filter.py
Identifies protein-coding genes from the processed protein FASTA file and retains only the corresponding gene information extracted from the GFF file.

Identify_neighboring_genes_and_calculate_intergenic_distancey.py
Identifies the upstream and downstream neighboring genes of each gene and retrieves their associated information and intergenic distances.

3. Calculate codon composition and relative synonymous codon usage (RSCU) values

ATGC_calculate.py
Computes genome-wide codon usage, including third-position base frequencies and GC content at the first, second, and third codon positions (GC1s, GC2s, GC3s).

RSCU_calculate.py
Calculates the relative synonymous codon usage (RSCU) values for 61 codons (exclude start and stop codons) in the genome.

4. Concatenate aligned single-copy orthologous genes

Single_copy_orthologous_MSA_concatenated.py
Concatenates multiple sequence alignments of individual single-copy orthologous genes into a single alignment matrix for phylogenetic tree construction。


