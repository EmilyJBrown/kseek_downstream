## take a list of all files with readID and kmer, and write a file with counts of readpairs with ends in each combo of kmers

import sys
import gzip
import subprocess as sp
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

filelistfh=sys.argv[1]
newfh=sys.argv[2]

filelist=open(filelistfh, 'r')
newfile=open(newfh, 'w')

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

nfiles=0
kmerdict={}
kcounts={}

for line in filelist:
	if line=='': break
	nfiles+=1
	line=line.strip()
	myfh=line.split()[0]
	myfile=gzip.open(myfh)
	readdict={}
## first, populate a dictionary with read:[kmers]
	for line in myfile:
		if line=='': break
		line=line.strip()
		split=line.split()
		readid=split[0]
		kmer=split[1]
		kmer=simplekmer(kmer)
		if readid in readdict: readdict[readid].append(kmer)
		if readid not in readdict: readdict[readid]=[kmer]
		if kmer in kcounts: kcounts[kmer]+=1
		if kmer not in kcounts: kcounts[kmer]=1
## now go through that dictionary and find reads with both ends IDed as kmers
	for read in readdict:
		if len(readdict[read])<2: continue
		if len(readdict[read])>1:
			k1=readdict[read][0]
			k2=readdict[read][1]
			if k1 in kmerdict:
				if k2 in kmerdict[k1]: kmerdict[k1][k2]+=1
				if k2 not in kmerdict[k1]: kmerdict[k1][k2]=1
			if k1 not in kmerdict: kmerdict[k1]={k2:1}

filelist.close()
print "Done reading in all the kmers!"

## now go throught kmerdict, discard kmers with counts less than 100/ sample on average, then make a new dict 
## with 0 for all kmer pairs, then loop through dictionary and add to those

finalks=[]

for kmer in kcounts:
	mycount=kcounts[kmer]
	if mycount<=100*nfiles: continue
	if mycount>100*nfiles: finalks.append(kmer)

matdict={}
for j in sorted(finalks):
	if j not in matdict: matdict[j]={}
	for i in sorted(finalks):
		if i not in matdict[j]: matdict[j][i]=0

for k1 in kmerdict:
	if k1 not in matdict: continue
	if k1 in matdict:
		for k2 in kmerdict[k1]:
			if k2 not in matdict: continue
			if k2 in matdict:
				matdict[k1][k2]+=kmerdict[k1][k2]
				matdict[k2][k1]+=kmerdict[k1][k2]

for kmer in sorted(finalks): newfile.write(str('\t')+str(kmer))
newfile.write('\n')

for j1 in sorted(finalks):
	newfile.write(j1)
	for j2 in sorted(finalks):
		if j2==j1: newfile.write('\t'+str(0))
		if j2!=j1:
			newfile.write('\t'+str(float(matdict[j1][j2])/float(kcounts[j1]+kcounts[j2])))
	newfile.write('\n')

newfile.close()



