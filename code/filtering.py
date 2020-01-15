# filtering.py

# Filtering the matrix 

class Filter(object):
    def process(self,data_list,kmer_idx,strain_idx):
        raise NotImplemented()

# Simplest possible filter, only kmers that appear in:
# at least min_presence times
# at most max_presence times
class SimpleFilter(Filter):
    def __init__(self,min_presence,max_presence):
        self.min_presence = min_presence
        self.max_presence = max_presence

    def process(self,data_list,kmer_idx,strain_idx):
        retain = set()
        reverse_kmer_idx = {v:k for k,v in kmer_idx.items()}

        for kmer,pos in kmer_idx.items():
            entries = list(filter(lambda x: x[0] == pos),data_list)
            n_entries = len(entries)
            if n_entries >= self.min_presence and n_entries <= self.max_presence:
                retain.add(kmer) # kmer position to retain

        
        # make a new strain index
        new_kmer_idx = {}
        new_kmer_pos = 0
        new_data_list = []
        for row,col,count in data_list:
                kmer = reverse_kmer_idx[row]
                if kmer in remain:
                    if not kmer in new_kmer_idx:
                        new_kmer_idx[kmer] = new_kmer_pos
                        new_kmer_pos += 1

                        new_row = new_kmer_idx[kmer]
                        new_data_list.append((new_row,col,count))