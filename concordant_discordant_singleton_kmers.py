## takes a list of files with readIDs and kmer identified from kseek from multiple samples
## writes a file, kmer_pairs_singletons.txt, with counts of each observed pair of kmer for each sample

import sys
import gzip
import subprocess as sp
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def revcomp(kmer):
	mydna=Seq(kmer, generic_dna)
	myrc=mydna.reverse_complement()
	return myrc

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

filelistfh=sys.argv[1]
filelist=open(filelistfh, 'r')

tempfile=open("tempfile.tmp", "w")

allsamples={}

kpairs={}

newfile=open("kmer_pairs_singletons.txt", "w")

for line in filelist:
	if line=='': break
	line=line.strip()
	myfile=gzip.open(line, 'r')
	sample=line.split("_")[0]
	if sample not in allsamples: allsamples[sample]=1
	readdict={}
	for aline in myfile:
		if aline=='': break
		aline=aline.strip()
		readid=aline.split()[0]
		kmer=aline.split()[1]
		kmer=simplekmer(kmer)
		if readid in readdict: readdict[readid].append(kmer)
		if readid not in readdict: readdict[readid]=[kmer]
## now reorganize them by kmer pairs, i.e. either k1:k1, k1:k2, or k1:none, and load into kpairs dictionary	
	for myread in readdict:
## first get the singletons
		if len(readdict[myread])==1: 
			myks=(readdict[myread][0],"none")
			if myks not in kpairs: kpairs[myks]={}
			if sample in kpairs[myks]: kpairs[myks][sample]+=1
			if sample not in kpairs[myks]: kpairs[myks][sample]=1
## then get all the reads where both ends have some kmer identified, and alphabetize them since we don't care which came first
		if len(readdict[myread])==2:
			myks=sorted(readdict[myread])
			myks=(myks[0], myks[1])
			if myks not in kpairs: kpairs[myks]={}
			if sample in kpairs[myks]: kpairs[myks][sample]+=1
			if sample not in kpairs[myks]: kpairs[myks][sample]=1
	myfile.close()

## now go through the kpairs dictionary and write to the new file!

newfile.write("Kmer1\tKmer2\t")
newfile.write("\t".join(sorted(allsamples.keys())))
newfile.write('\n')

allpairs=kpairs.keys()

for myp in sorted(allpairs):
	tempfile.write(str(myp)+'\n')

for apair in sorted(allpairs):
	newfile.write(str(apair[0]   )+'\t'+str(apair[1]))
	for s in sorted(allsamples.keys()):
		if s not in kpairs[apair]: newfile.write("\t0")
		if s in kpairs[apair]: newfile.write("\t%i" %(kpairs[apair][s]))
	newfile.write('\n')

filelist.close()
newfile.close()



