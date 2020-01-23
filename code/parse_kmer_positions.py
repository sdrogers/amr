import argparse
import json
import os


def filename_to_strain(filename):
    base = os.path.basename(filename)
    strain = '_'.join(base.split('_')[0:2])
    return strain


def process_file(filename):
    with open(filename, 'rb') as f:
        raw_data = json.load(f)

    # enable file name lookup by index
    # note - we're not verifying that the indices are unique!
    file_lookup = dict([x[::-1] for x in raw_data['file_idx'].items()])

    for file_idx, contig, pos, kmer_id_strand in raw_data['position_list']:
        kmer_file_name = file_lookup[file_idx]
        kmer_strain = filename_to_strain(kmer_file_name)
        kmer_id = abs(kmer_id_strand)

        output_string = '\t'.join([str(x) for x in (kmer_strain, contig, kmer_id, pos, pos + 31)])
        print(output_string)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    args = parser.parse_args()

    filename = args.file
    process_file(filename)

