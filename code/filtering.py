# filtering.py

# Filtering the matrix 

class Filter(object):
    def process(self,data_list,kmer_idx,strain_idx):
        raise NotImplemented()

# Simplest possible filter, only kmers that appear in:
# at least min_presence times
# at most max_presence times
class SimpleFilter(Filter):
    def __init__(self,min_presence = 5,max_presence = 1000):
        self.min_presence = min_presence
        self.max_presence = max_presence

    def process(self,data_list,kmer_idx,strain_idx):
        retain = set()
        reverse_kmer_idx = {v:k for k,v in kmer_idx.items()}

        kmer_counts = {k:0 for k in kmer_idx}
        for kmer_pos,strain_pos,count in data_list:
            kmer = reverse_kmer_idx[kmer_pos]
            kmer_counts[kmer] += 1        


        for kmer,count in kmer_counts.items():
            if count >= self.min_presence and count <= self.max_presence:
                retain.add(kmer)
        
        # make a new kmer index
        new_kmer_idx = {}
        new_kmer_pos = 0
        new_data_list = []
        for row,col,count in data_list:
                kmer = reverse_kmer_idx[row]
                if kmer in retain:
                    if not kmer in new_kmer_idx:
                        new_kmer_idx[kmer] = new_kmer_pos
                        new_kmer_pos += 1

                        new_row = new_kmer_idx[kmer]
                        new_data_list.append((new_row,col,count))
        
        return new_data_list,new_kmer_idx,strain_idx