## takes a file with a list of bam alignments, returns a file with the distance between mate pairs for uniquely mapped sequences
## mate pairs must also map to the same chromosome

import sys
import gzip
import subprocess as sp
import os

filelistfh=sys.argv[1]
filelist=open(filelistfh, 'r')

for line in filelist:
	if line=='': break
	readdict={}
	bamfh=line.strip()
	proc1=sp.Popen(['samtools','view',bamfh], stdout=sp.PIPE)
	output1=proc1.stdout
	newfh=bamfh.split('/')[-1]
	newfh=newfh.split("_")[0]
	newfh=str(newfh)+str("_distance_between_matepairs.txt.gz")
	newfile=gzip.open(newfh, 'wb')
	for aline in output1:
		if aline=='': break
		aline=aline.strip()
		split=aline.split()
		readid=split[0]
		chrom=split[2]
		pos=int(split[3])
		uniq=split[12]
		if uniq.startswith("XS:i"): continue
		if readid in readdict: 
			if chrom in readdict[readid]:
				newfile.write(str(chrom)+'\t'+str(abs(pos-readdict[readid].pop(chrom)))+'\n')
			if chrom not in readdict[readid]: continue
		if readid not in readdict: readdict[readid]={chrom:pos}
	print(str("Got all the distances for sample ")+str(newfh))
	newfile.close()

filelist.close()

