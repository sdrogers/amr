#!/usr/bin/env python
# coding: utf-8

from statsmodels.stats import multitest
import argparse
import scipy.stats


def hg_test(tagged_objects, sampled_objects, object_count):
    tagged_objects = set(tagged_objects)
    sampled_objects = set(sampled_objects)

    intersection = len(tagged_objects.intersection(sampled_objects))
    total_tagged = len(tagged_objects)
    sample_size = len(sampled_objects)

    M = object_count
    n = total_tagged
    N = sample_size
    k = intersection

    return scipy.stats.hypergeom.sf(k, M, n, N, 1)


def process_family_file(filename, associations, contig_list, adjust=True):
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

    id_list = []
    pvalues = []
    for kmer_id, associated_bgcs in associations.items():
        # strain_id, contig_id, kmer_id = kmer_identifier
        kmer_bgc_ids = ["{}.cluster{:03}".format(cont, int(x)) for strain, cont, x in associated_bgcs]

        for gcf_id, gcf_bgc_ids in gcfs.items():
            prob = hg_test(kmer_bgc_ids, gcf_bgc_ids, total_bgc_count)

            id_list.append((kmer_id, gcf_id))
            pvalues.append(prob)

    if adjust:
        res = multitest.multipletests(pvalues)
        reject = res[0]
        pvalues = res[1]

    for (kmer_id, gcf_id), prob in zip(id_list, pvalues):
        yield kmer_id, gcf_id, prob


def load_bgc_kmer_associations(bgc_kmers_file):
    associated_bgcs_by_kmer = {}
    with open(bgc_kmers_file, 'r') as f:
        for l in f:
            strain_id, contig_id, bgc_id, bgc_start, bgc_end, kmer_id, kmer_start, kmer_end = l.strip().split()
            if kmer_id not in associated_bgcs_by_kmer:
                associated_bgcs_by_kmer[kmer_id] = set([])
            associated_bgcs_by_kmer[kmer_id].add((strain_id, contig_id, bgc_id))

    # If the same set of BGCs gets tagged by several kmers, only choose one of them
    # This should probably be done intelligently!

    kmer_assoc = []
    bgc_assoc = []
    for kmer_id, bgc_id_set in associated_bgcs_by_kmer.items():
        bgc_ids = list(bgc_id_set)
        bgc_ids.sort()
        if bgc_ids not in bgc_assoc:
            bgc_assoc.append(bgc_ids)
            kmer_assoc.append(set([]))
        idx = bgc_assoc.index(bgc_ids)
        kmer_assoc[idx].add(kmer_id)

    # Do the back translation, selecting a representative kmer where they tag identical things

    associated_bgcs_by_kmer_filtered = {}
    for bgc_ids, kmer_ids in zip(bgc_assoc, kmer_assoc):
        kmer_ids = list(kmer_ids)
        kmer_ids.sort()
        kmer_id = kmer_ids[0]
        associated_bgcs_by_kmer_filtered[kmer_id] = bgc_ids

    return associated_bgcs_by_kmer_filtered


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

    associated_bgcs_by_kmer = load_bgc_kmer_associations(bgc_kmers_file)

    print("# kmer_id\tgcf_id\tp_val")
    for kmer_id, gcf_id, prob in process_family_file(gcf_file, associated_bgcs_by_kmer, contig_list):
        if prob <= threshold:
            print("{}\t{}\t{}".format(kmer_id, gcf_id, prob))






