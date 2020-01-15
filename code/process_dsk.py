# process_dsk.py
import csv

class DSKProcess(object):
    def __init__(self,filename_list):
        self.filename_list = filename_list
    
    def process(self):
        kmer_idx = {}
        kmer_pos = 0
        stain_idx = {}
        strain_pos = 0
        for i,filename in enumerate(filename_list):
            print(filename, "({} of {})".format(i,len(filename_list)))
            strain_name = Path(filename).stem.split('_')[1] # FIX ME
            with open(filename,'r') as f:
                reader = csv.reader(f)
                heads = next(reader) # check it has a heading row
                for line in reader:
                    kmer,count = line
                    if not kmer in kmer_idx:
                        kmer_idx[kmer] = kmer_pos
                        kmer_pos += 1
                    
                    strain_idx[strain_name] = strain_pos
                    strain_pos += 1

                    row = kmer_idx[kmer]
                    col = strain_idx[strain_name]

                    data_list.append((row,col,count))
        
        return data_list,kmer_idx,strain_idx

