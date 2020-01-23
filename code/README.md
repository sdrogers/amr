python parse_kmer_positions.py <example_data>/kmer_positions_pyraz.json > <example_data>/strain_kmer_positions.tsv
python collect_antismash_bgcs.py <antismash_folder> > <example_data>/bgc_positions.tsv
python compare_positions.py <example_data>/bgc_positions.tsv <example_data>/strain_kmer_positions.tsv -d <margin> > cooccurring_bgcs_strains.tsv

