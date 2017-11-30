## takes a bamfile and a gzip-ed fastq file, writes a new fastq file with reads that are NOT in the bam file

import sys
import gzip
import subprocess as sp
import os

bamfh=sys.argv[1]
fastqfh=sys.argv[2]
newfh=sys.argv[3]

fastqfile=gzip.open(fastqfh, 'r')
newfile=gzip.open(newfh, 'wb')

mapdict={}

proc1=sp.Popen(['samtools','view',bamfh], stdout=sp.PIPE)
output1=proc1.stdout
for line in output1:
	if line=='': break
	line=line.strip()
	split=line.split()
	readID=split[0]
	if readID not in mapdict: mapdict[readID]=1

fastqcount=0
unmapped=0
write="no"
myid=0
for line in fastqfile:
	fastqcount+=1
	if fastqcount%4==1: 
		myid=line.split()[0]
		if myid in mapdict: 
			write="no"
		if myid not in mapdict: 
			write="yes"
			unmapped+=1
	if fastqcount%4==2:
		if write=="yes":
			newfile.write(str(">")+str(myid)+'\n')
			newfile.write(str(line))
		if write=="no": 
			continue
	else:
		continue
fastqfile.close()
newfile.close()
