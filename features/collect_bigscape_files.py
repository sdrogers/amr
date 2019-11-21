#!/usr/bin/python

import os


def collect_bigscape_files(bigscape_path):
    """
    Collect the results of a bigscape run (curently hardcoded to cutoff .30)
    output them to stdout in csv format
    columns are:
        family_id
        bgc_id
        strain
    """
    families = {}
    
    for prod_type in os.listdir(bigscape_path):
        if prod_type.endswith('.tsv'):
            continue
        filename = '{}/{}/{}_clustering_c0.30.tsv'.format(bigscape_path, prod_type, prod_type)
        with open(filename) as f:
            for l in f.readlines():
                if l.startswith('#'):
                    continue
                bgc, family = l.strip().split()
                family_name = prod_type + '_' + family
                bgc_id = bgc.split('.')[-1]
                bgc_idx = str(int(bgc_id[-3:]))
                strain_id = '.'.join(bgc.split('.')[:-1])
    
                if family_name in families:
                    families[family_name].append((strain_id, bgc_idx))
                else:
                    families[family_name] = [(strain_id, bgc_idx)]
    
    print('#family_id,bgc_id,strain')
    for fam_id, members in families.items():
        for strain, bgc_id in members:
            print('{},{},{}'.format(fam_id, bgc_id, strain))


if __name__ == '__main__':
    bigscape_path = "/srv/data/tuberculosis/bigscape/20190621-bigscape-auto/network_files/2019-09-03_17-16-46_hybrids_auto"
    collect_bigscape_files(bigscape_path)
