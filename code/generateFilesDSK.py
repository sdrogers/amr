import sys
import os
import subprocess
import glob

class KmerFiles (object):
	#might want to change some of the other dsk parameters at some point ?
	def __init__(self,directoryGenomes,dsk_path,kmer_abundnace_min,kmer_size):
		self.directoryGenomes = directoryGenomes
		self.dsk_path = dsk_path
		self.kmer_abundnace_min = kmer_abundnace_min
		self.kmer_size = kmer_size

	def process(self):
		genomeFiles = glob.glob(self.directoryGenomes+"/"+ "*.fna")
		print (len(genomeFiles))
		for f in genomeFiles:
			subprocess.check_output(f'{self.dsk_path}/bin/dsk -file {f} -abundance-min {self.kmer_abundnace_min} -kmer-size {self.kmer_size} -out-dir {self.directoryGenomes}/kmersh5',shell=True)

		for f in genomeFiles:
			fn,ext = os.path.splitext(os.path.basename(f))
			subprocess.check_output(f'{self.dsk_path}/bin/dsk2ascii -file {self.directoryGenomes}/kmersh5/{fn}.h5 -out {self.directoryGenomes}/kmersh5/{fn}.txt', shell=True)

		return glob.glob(self.directoryGenomes+"/kmersh5/"+"*.txt")



if __name__ == '__main__':
	ks = KmerFiles("/Users/alexpancheva/Documents/dskTESTFiles",
		"/Users/alexpancheva/Documents/PATRIC/dsk-v2.3.0-bin-Darwin",1,31)
	print (ks.process())





