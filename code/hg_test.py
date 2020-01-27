#!/usr/bin/env python
# coding: utf-8

import argparse
import scipy.stats


def hg_test(tagged_strains, sampled_strains, strain_count):
    tagged_strains = set(tagged_strains)
    sampled_strains = set(sampled_strains)

    intersection = len(tagged_strains.intersection(sampled_strains))
    total_tagged = len(tagged_strains)
    sample_size = len(sampled_strains)

    M = strain_count
    n = total_tagged
    N = sample_size
    k = intersection

    return scipy.stats.hypergeom.sf(k, M, n, N, 1)


def process_family_file(filename, associations, contig_list):
    gcfs = {}
    all_bgcs = set([])

    contig_list = [x.split('.')[0] for x in contig_list]
    with open(filename, 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue
            bgc_id, gcf_id = l.strip().split()
            if bgc_id.split('.')[0] not in contig_list:
                continue
            if gcf_id not in gcfs:
                gcfs[gcf_id] = set([])
            gcfs[gcf_id].add(bgc_id)
            all_bgcs.add(bgc_id)

    total_bgc_count = len(all_bgcs)

    for kmer_identifier, associated_bgcs in associations.items():
        strain_id, contig_id, kmer_id = kmer_identifier
        kmer_bgc_ids = ["{}.cluster{:03}".format(contig_id, int(x)) for x in associated_bgcs]

        for gcf_id, gcf_bgc_ids in gcfs.items():
            prob = hg_test(kmer_bgc_ids, gcf_bgc_ids, total_bgc_count)

            yield kmer_id, gcf_id, prob



# kmers_file = "/home/grimur/git/amr/example_data/strain_kmer_positions.tsv"
# bgc_file = "/home/grimur/git/amr/example_data/bgc_positions.tsv"
# bgc_kmers_file = "/home/grimur/git/amr/example_data/bgc_kmers.tsv"
# bigscape_data_folder = "/home/grimur/uist-data/tuberculosis/bigscape/20190621-bigscape-auto/network_files/2019-09-03_17-16-46_hybrids_auto/"
# nrps_gcfs_filename = bigscape_data_folder + 'NRPS/NRPS_clustering_c0.50.tsv'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate hypergeometric correlation btw. GCF and k-mer')
    parser.add_argument('-g', dest='gcf', help='gcf file', required=True)
    parser.add_argument('-k', dest='kmer', help='kmer file', required=True)
    parser.add_argument('-c', dest='bgc_kmer', help='kmer-bgc correlationfile', required=True)
    parser.add_argument('-t', dest='thresh', help='p value threshold', type=float, default=1.0)
    args = parser.parse_args()

    gcf_file = args.gcf
    bgc_kmers_file = args.bgc_kmer
    kmer_file = args.kmer
    threshold = args.thresh

    strain_contigs = set([])
    with open(kmer_file, 'r') as f:
        for l in f:
            strain_id, contig_id, amr_id, start, end = l.strip().split()
            strain_contigs.add((strain_id, contig_id))

    contig_list = [x[1] for x in strain_contigs]

    associated_bgcs_by_kmer = {}
    with open(bgc_kmers_file, 'r') as f:
        for l in f:
            strain_id, contig_id, bgc_id, bgc_start, bgc_end, kmer_id, kmer_start, kmer_end = l.strip().split()
            if (strain_id, contig_id, kmer_id) not in associated_bgcs_by_kmer:
                associated_bgcs_by_kmer[(strain_id, contig_id, kmer_id)] = set([])
            associated_bgcs_by_kmer[(strain_id, contig_id, kmer_id)].add(bgc_id)

    print("# kmer_id\tgcf_id\tp_val")
    for kmer_id, gcf_id, prob in process_family_file(gcf_file, associated_bgcs_by_kmer, contig_list):
        if prob <= threshold:
            print("{}\t{}\t{}".format(kmer_id, gcf_id, prob))






