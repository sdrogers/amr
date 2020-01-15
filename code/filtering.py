# filtering.py

# Filtering the matrix 

class Filter(object):
    def process(self,input_matrix,kmer_list,strain_list):
        raise NotImplemented()

class SimpleFilter(Filter):
    def process(self,input_matrix,kmer_list,strain_list):
        