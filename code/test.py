#test.py

from process_dsk import DSKProcess

d = DSKProcess(['/Users/simon/git/amr/example_data/dsk/dskTxtExample.txt'])
data_list,kmer_idx,strain_idx = d.process()
print(len(data_list))

from filtering import SimpleFilter
s = SimpleFilter(min_presence = 1)
new_data_list,new_kmer_idx,stain_idx = s.process(data_list,kmer_idx,strain_idx)

print(len(new_data_list))