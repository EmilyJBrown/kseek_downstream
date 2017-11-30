## takes a file with sets of files, with the gzip-ed readIDs first, then the two paired read gzip-ed fastq file
## write two new gzip-ed fastq files, forward and reverse, with the reads indicated in the first file

import sys
import gzip

filelistfh=sys.argv[1]
filelist=open(filelistfh, 'r')

for line in filelist:
	readsdict={}
	if line=='': break
	line=line.strip()
	split=line.split()
	readidsfh=split[0]
	readids=gzip.open(readidsfh, 'r')
	for aline in readids:
		if aline=='': break
		aline=aline.strip()
		asplit=aline.split()
		myid=asplit[0]
		if myid not in readsdict: readsdict[myid]=1
	readids.close()
	print "Got my read IDs from file"
	print readidsfh
	freadsfh=split[1]
	freads=gzip.open(freadsfh, 'r')
	readsbase=freadsfh.split('.')[0]
	readsbase=readsbase.split('/')[-1]
	newfreadsfh=str(readsbase)+str('.IDed_kseek.fq.gz')
	newfreads=gzip.open(newfreadsfh, 'wb')
	linecount=0
	write="no"
	for bline in freads:
		if bline=='': break
		bline=bline.strip()
		linecount+=1
		if linecount%4==1:
			blineid=bline.split()[0]
			if blineid in readsdict:
				write="yes"
				newfreads.write(str(bline)+'\n')
			if blineid not in readsdict:
				write="no"
				continue
		if linecount%4!=1 and write=="yes": newfreads.write(str(bline)+'\n')
		if linecount%4!=1 and write=="no": continue
	newfreads.close()
	freads.close()
	print "Wrote all the forward reads!"
	rreadsfh=split[2]
	rreads=gzip.open(rreadsfh, 'r')
	readsbase=rreadsfh.split('.')[0]
	readsbase=readsbase.split('/')[-1]
	newrreadsfh=str(readsbase)+str('.IDed_kseek.fq.gz')
	newrreads=gzip.open(newrreadsfh, 'wb')
	linecount=0
	write="no"
	for cline in rreads:
		if cline=='': break
		cline=cline.strip()
		linecount+=1
		if linecount%4==1:
			clineid=cline.split()[0]
			if clineid in readsdict:
				write="yes"
				newrreads.write(str(cline)+'\n')
			if clineid not in readsdict:
				write="no"
				continue
		if linecount%4!=1 and write=="yes": newrreads.write(str(cline)+'\n')
		if linecount%4!=1 and write=="no": continue
	newrreads.close()
	rreads.close()
	print "Wrote all the reverse reads! Moving on to next file."
filelist.close()
