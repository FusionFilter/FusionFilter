# Instructions for computing approximate paralog clusters and simpler blast match clusters for annotating suspicious fusion calls.


## Sort blast pairs

    gunzip -c  blast_pairs.gene_syms.outfmt6.gz | sort -k4,4g -k3,3gr > blastn.outfmt6.grouped.geneSym.sorted


## Pull just the top hits:

    get_top_blast_pairs.pl blastn.outfmt6.grouped.geneSym.sorted > blastn.outfmt6.grouped.geneSym.sorted.top


## perform paralog-level clustering:

   outfmt6_add_percent_match_length.group_segments.to_Markov_Clustering.pl --outfmt6_grouped blastn.outfmt6.grouped.geneSym.sorted.top --min_pct_len 1 --min_per_id 90 --inflation_factor 5


## perform simple 'nucleotide blast clusters'

    outfmt6_add_percent_match_length.group_segments.to_Markov_Clustering.pl --outfmt6_grouped blastn.outfmt6.grouped.geneSym.sorted.top --min_pct_len 1 --min_per_id 60 --inflation_factor 5

