import argparse


def parse_file(filename):
    with open(filename, 'r') as f:
        for l in f.readlines():
            strain_id, contig_id, object_id, start, stop = l.strip().split()
            start = int(start)
            stop = int(stop)
            yield strain_id, contig_id, object_id, start, stop


def parse_strain_from_file(strain_id, filename):
    for l in parse_file(filename):
        if l[0] == strain_id:
            yield l


def file_len(filename):
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def load_reference_data(filename):
    """
    Reference data is structured by strain id and contig id in nested
    dictionaries with positions as a list of tuples:

    reference = {
        strain_1: {
            contig_1: [
                (start, stop),
                (start, stop),
                ...
            ],
            contig_2: [
                (start, stop),
                ...
            ]
        },
        strain_1: {
            ...
        },
        ...
    }
    """
    reference = {}
    for strain_id, contig_id, object_id, start, stop in parse_file(filename):
        if strain_id not in reference:
            reference[strain_id] = {}
        if contig_id not in reference[strain_id]:
            reference[strain_id][contig_id] = []

        reference[strain_id][contig_id].append((start, stop, object_id))

    # sort the positions
    for strain_id in reference.keys():
        for contig_id in reference[strain_id].keys():
            reference[strain_id][contig_id].sort(key=lambda x: x[0])

    return reference


def check_overlap_seq(start, stop, pos_list, margin):
    for start_pos, stop_pos, obj_id in pos_list:
        if stop > start_pos - margin and start < stop_pos + margin:
            yield start_pos, stop_pos, obj_id

        # Looking at stuff beyond the end of the section
        if start_pos > stop + margin:
            break


def compare_files(file_1, file_2, margin):
    # Use the smaller file as reference and loop through the larger one
    if file_len(file_1) < file_len(file_2):
        # Keep output order same as input
        reverse_order = True
        base_file = file_1
        check_file = file_2
    else:
        reverse_order = False
        base_file = file_2
        check_file = file_1

    reference_data = load_reference_data(base_file)

    for strain_id, contig_id, object_id, start, stop in parse_file(check_file):
        if strain_id in reference_data:
            if contig_id in reference_data[strain_id]:
                for start_pos, stop_pos, obj_id in check_overlap_seq(start, stop, reference_data[strain_id][contig_id], margin):
                    if reverse_order:
                        output_fields = (strain_id, contig_id, obj_id, start_pos, stop_pos, object_id, start, stop)
                    else:
                        output_fields = (strain_id, contig_id, object_id, start, stop, obj_id, start_pos, stop_pos)
                    output = '\t'.join([str(x) for x in output_fields])
                    print(output)


def test_check_overlap_seq():
    pos_list = [(10, 20, 'a'), (30, 40, 'b'), (50, 60, 'c'), (70, 80, 'd')]

    for start_pos, stop_pos, obj_id in check_overlap_seq(25, 26, pos_list):
        print(obj_id)
    print()
    for start_pos, stop_pos, obj_id in check_overlap_seq(45, 46, pos_list):
        print(obj_id)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('f1')
    parser.add_argument('f2')
    parser.add_argument('-d', dest='margin', help='distance btw. BGC and k-mer', type=int, required=False, default=10)
    args = parser.parse_args()

    file1 = args.f1
    file2 = args.f2
    margin = args.margin

    compare_files(file1, file2, margin)
