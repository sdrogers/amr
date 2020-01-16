# process_dsk.py
import csv
from pathlib import Path



class DSKProcess(object):
    def __init__(self,filename_list):
        self.filename_list = filename_list
    
    def process(self):
        kmer_idx = {}
        kmer_pos = 0
        strain_idx = {}
        strain_pos = 0
        data_list = []
        for i,filename in enumerate(self.filename_list):
            print(filename, "({} of {})".format(i,len(self.filename_list)))
            # Just use the filename at this point....

            # try: # fix for the test data - what format will filenames be in??
            #     strain_name = Path(filename).stem.split('_')[1]
            # except:
            #     strain_name = filename
            strain_name = filename
            
            with open(filename,'r') as f:
                strain_idx[strain_name] = strain_pos
                strain_pos += 1
                reader = csv.reader(f,delimiter = ' ')
                for line in reader:
                    kmer,count = line
                    if not kmer in kmer_idx:
                        kmer_idx[kmer] = kmer_pos
                        kmer_pos += 1
                    
                    row = kmer_idx[kmer]
                    col = strain_idx[strain_name]

                    data_list.append((row,col,int(count)))
        
        return data_list,kmer_idx,strain_idx


if __name__ == '__main__':
    d = DSKProcess(['/Users/simon/git/amr/example_data/dsk/dskTxtExample.txt'])
    data_list,kmer_idx,strain_idx = d.process()
    print(data_list)