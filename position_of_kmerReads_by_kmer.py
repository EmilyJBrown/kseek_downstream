## take a filelist with bam files and gziped readIDs and kmer, 
## write a new file with kmers mapped at each position in the genome

import sys
import gzip
import subprocess as sp
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

myfh=sys.argv[1]
newfh=sys.argv[2]
myfile=open(myfh, 'r')
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


posdict={}

for line in myfile:
	if line=='': break
	line=line.strip()
	bamfh=line.split()[0]
	idsfh=line.split()[1]
	idsfile=gzip.open(idsfh)
	kmerdict={}
	for line in idsfile:
		if line=='': break
		line=line.strip()
		readid=line.split()[0]
		readid=readid.split("@")[1]
		kmer=line.split()[1]
		alphakmer=simplekmer(kmer)
		if readid not in kmerdict: kmerdict[readid]=alphakmer
	idsfile.close()
	proc1=sp.Popen(['samtools','view',bamfh], stdout=sp.PIPE)
	output1=proc1.stdout
	for line in output1:
		if line=='': break
		line=line.strip()
		split=line.split()
		myid=split[0]
		if myid not in kmerdict: continue
		mykmer=kmerdict[myid]
		chrom=split[2]
		pos=int(split[3])
		if mykmer in posdict:
			if chrom in posdict[mykmer]:posdict[mykmer][chrom].append(pos)
			if chrom not in posdict[mykmer]: posdict[mykmer][chrom]=[pos]
		if mykmer not in posdict: posdict[mykmer]={chrom:[pos]}
	print str("Got the positions for sample ")+str(bamfh)

print "Done getting all the positions for all the kmers.  Writing to new file."

newfile.write(str("kmer")+'\t'+str("chrom")+'\t'+str("pos")+'\n')
allkmers=posdict.keys()
for k in sorted(allkmers):
	mychroms=posdict[k]
	for c in sorted(mychroms):
		mypos=posdict[k][c]
		for p in sorted(mypos):
			newfile.write(str(k)+'\t'+str(c)+'\t'+str(p)+'\n')
newfile.close()

