## take a list of kmers that you've identified as high copy number from k-seek
## writes a file with ~100bp of tandem repeat of all permutations of those kmers in fasta format

import sys
import gzip
import subprocess as sp
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def revcomp(kmer):
	mydna=Seq(kmer, generic_dna)
	myrc=mydna.reverse_complement()
	return myrc

def allkmers(kmer):
	subks=[]
	doublek=str(kmer)+str(kmer)
	myrc=revcomp(kmer)
	doublerc=str(myrc)+str(myrc)
	mylen=len(kmer)
	for i in range(0,mylen):
		mysub1=doublek[i:i+mylen]
		subks.append(mysub1)
		mysub2=doublerc[i:i+mylen]
		subks.append(mysub2)
	mysubks=sorted(subks)
	return mysubks

def simplekmer(kmer):
	subks=[]
	doublek=str(kmer)+str(kmer)
	myrc=revcomp(kmer)
	doublerc=str(myrc)+str(myrc)
	mylen=len(kmer)
	mysubs=[]
	for i in range(0,mylen):
		mysub1=doublek[i:i+mylen]
		subks.append(mysub1)
		mysub2=doublerc[i:i+mylen]
		subks.append(mysub2)
	mysubks=sorted(subks)
	firstk=mysubks[0]
	return firstk

kfh=sys.argv[1]
newfh=sys.argv[2]

kfile=open(kfh, 'r')
newfile=open(newfh, 'w')

for line in kfile:
	if line=='': break
	line=line.strip()
	if line.startswith("lines"): continue
	myk=line.split()[0]
	myk=myk.split('/')[0]
	mylen=len(myk)
	if mylen<2: continue
	myrep=int(100/mylen)
	newfile.write(str(">")+str(myk)+'\n')
	newfile.write(str(myrep*myk)+'\n')

kfile.close()
newfile.close()