## run this from a directory with the .rep files output from k-seek
## it will get read IDs for k-seek reads for each sample, merging split k-seek outputs 
## it will ignore reads that have a 1-mer identified

import sys
import gzip
import os

allfiles=os.listdir('./')
myfiles=[]
for file in allfiles:
	if file.endswith("rep"): myfiles.append(file)
print myfiles

readsdict={}

for afile in myfiles:
	print afile
	myid=afile.split('_')[0]
	print myid
	if myid not in readsdict: readsdict[myid]={}
	myfile=open(afile, 'r')
	count=0
	readid=0
	for line in myfile:
		count+=1
		if count%5==1:
			readid=line.split()[0]
		if count%5==0:
			kmer=line.split("=")[0]
			if len(kmer)>1: 
				if readid in readsdict[myid]: readsdict[myid][readid].append(kmer)
				if readid not in readsdict[myid]: readsdict[myid][readid]=[kmer]
			if len(kmer)==1: continue
	myfile.close()

for lib in readsdict:
	mynewfh=str(lib)+str("_readIDs_from_kseek.txt.gz")
	mynewfile=gzip.open(mynewfh, 'wb')
	for readid in readsdict[lib]:
		for kmer in readsdict[lib][readid]:
			mynewfile.write(str(readid)+'\t'+str(kmer)+'\n')
	mynewfile.close()