#find_kmer_test.py

from find_kmers import FindKmers
import glob
import os
import random
import string
import pandas as pd
from pathlib import Path



def randomString(stringLength=10):
    """Generate a random string of fixed length """
    alphabet = ['A','C','G','T']
    return ''.join(random.choice(alphabet) for i in range(stringLength))

def randomSequences (n, stringLength = 10):
    """Generate a list of random sequences"""
    seqs = []
    for s in range(n):
        seqs.append(randomString(stringLength))
    return seqs

#kmer_list = randomSequences (10 ,stringLength = 7) 
kmer_list =[]
kmer_df = pd.read_csv ('prediction.csv')
for i,row in kmer_df.iterrows():
    if row.score>0:
        kmer_list.append(row.kmer)
print(len(kmer_list))

fasta_folder = '../data/fasta/'
all_filenames = glob.glob(os.path.join(fasta_folder,'*.fna'))    
ab = 'pyrazinamide'
mdata = pd.read_csv(f'../TBmetadata/{ab}.csv')
TB = (list(mdata[mdata.columns[0]]))
TB = ['000'+str(b) for b in TB]  # correct for - somewhere strain id was coverted to int ( removing 0s)
file_list = [f for f in all_filenames if Path(f).stem.split('_')[1] in TB]

print(f"number of kmers {len(kmer_list)}, number of files {len(file_list)}")
    
    
k = FindKmers(file_list[:5],kmer_list[:10])
position_list, kmer_idx, file_idx = k.process()
    
    
print(f"{len(position_list)} hits for: number of kmers {len(kmer_idx)}, number of files {len(file_idx)}")

