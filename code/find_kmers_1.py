# find_kmers.py
from Bio import SeqIO

class FindKmers(object):
    def __init__ (self, filename_list, kmer_list):
        self.filename_list = filename_list
        kmer_ids  = {k:i+1 for (i,k) in enumerate(kmer_list)}
        self.kmer_index = kmer_ids
                 
    
    def compliment(self,dna):
        """Return the reverse compliment of a sequence """
        bp = { 'A':'T', 'T':'A', 'C':'G','G':'C'}
        revdna =''    
        for i in range(len(dna)-1,-1,-1):
            revdna += bp[dna[i]]
        return (revdna)

    def find_all(self, a_str, sub):
        start = 0
        while True:
            start = a_str.find(sub, start)
            if start == -1: return
            yield start
            start += len(sub) # use start += 1 to find overlapping matches             
        
    def process(self):
        
        data_list = []        
        file_index = {f:i for (i,f) in enumerate(self.filename_list)}
                 
        for filename,file_pos in (file_index.items()):
            n_records = 0
             
            with open(filename, 'rt') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    contig_id =record.description.split()[0]
                    n_records += 1
                    seq =record.seq
                    print('****',file_pos, filename,contig_id, len(seq),n_records)
                
                    for kmer, kmer_id in self.kmer_index.items():
                        for pos in (list(self.find_all(seq,kmer))):   
                            data_list.append((file_pos,contig_id,pos,kmer_id))
                        for pos in (list(self.find_all(seq,self.compliment(kmer)))):
                            data_list.append((file_pos,contig_id,pos,-kmer_id)) 
        
        return data_list, self.kmer_index ,file_index
                                  
    
        
if __name__ == '__main__':
    k = FindKmers(['../data/fasta/GCF_000655135.1_Myco_tube_XTB13-183_V1_genomic.txt'],{'ACTCGGGCG': 0, 'AGCCCACTA': 1})
    position_list = k.process()
    print(position_list[:10])      
        
        
####