#test.py

from process_dsk import DSKProcess
import glob
import os

dsk_folder = '/srv/data/amr/sample_data/strain_kmers'
dsk_file_list = glob.glob(os.path.join(dsk_folder,'*.txt'))


d = DSKProcess(dsk_file_list)
data_list,kmer_idx,strain_idx = d.process()
print("n_kmers {}, n_strains {}".format(len(kmer_idx),len(strain_idx)))

from filtering import SimpleFilter
s = SimpleFilter(min_presence = round(0.05*len(kmer_idx)),max_presence = round(0.95*len(strain_idx)))
new_data_list,new_kmer_idx,stain_idx = s.process(data_list,kmer_idx,strain_idx)

from scipy.stats import coo_matrix

i,j,k = zip(*new_data_list)
coo = coo_matrix((k,(i,j)),shape = (len(new_kmer_idx),len(strain_idx)))
