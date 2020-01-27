# Parse the kmer json file to tsv format
python parse_kmer_positions.py <example_data>/kmer_positions_pyraz.json > <example_data>/strain_kmer_positions.tsv

# Parse the BGC positions from the antiSMASH result folder
python collect_antismash_bgcs.py <antismash_folder> > <example_data>/bgc_positions.tsv

# Find k-mers that occur in close proximity to BGCs
python compare_positions.py <example_data>/bgc_positions.tsv <example_data>/strain_kmer_positions.tsv -d <margin> > bgc_kmers.tsv

# Calculate probability that each GCF/k-mer pair has the actual or greater overlap in BGCs
python hg_test.py -g ../example_data/NRPS_clustering_c0.50.tsv -k ../example_data/strain_kmer_positions.tsv -c ../example_data/bgc_kmers.tsv -t 0.5
