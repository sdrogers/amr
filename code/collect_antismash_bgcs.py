import argparse
import os
from Bio import SeqIO


def collect_antismash_bgcs(path):
    """
    Parse antiSMASH output dir for BGCs.
    Print results to stout
    Results are
        strain_id
        contig_id
        cluster_id
        start
        stop
    """
    for strain_dir in os.listdir(path):
        antismash_output_path = os.path.join(path, strain_dir)
        for filename in os.listdir(antismash_output_path):
            if filename.endswith('.gbk') and 'final' in filename:
                file_path = os.path.join(antismash_output_path, filename)
                with open(file_path, 'r') as f:
                    for gb_record in SeqIO.parse(f, "genbank"):
                        strain_id = '_'.join(strain_dir.split('_')[:2])
                        contig_id = gb_record.id
                        for feature in gb_record.features:
                            if feature.type == 'cluster':
                                start = feature.location.start
                                end = feature.location.end
                                cluster_id = feature.qualifiers['note'][0].split()[-1]
                                output_list = [strain_id,contig_id, cluster_id, str(start), str(end)]
                                print('\t'.join(output_list))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    args = parser.parse_args()

    path = args.path
    collect_antismash_bgcs(path)
