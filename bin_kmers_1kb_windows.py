import sys

kmerfh=sys.argv[1]
newfh=sys.argv[2]

kmerfile=open(kmerfh, 'r')

newfile=open(newfh, 'w')

cdict={}
'''for line in chromfile:
	if line=='': break
	line=line.strip()
	split=line.split()
	chrom=split[0]
	clen=int(split[1])
	if chrom not in cdict: cdict[chrom]={}
	for i in range(0,clen):
		if i%1000==0:
			mywin=(i,i+1000)
			if mywin not in cdict[chrom]: cdict[chrom][mywin]={}
chromfile.close()
print "Got all my windows. Counting where the kmers are."
'''
for line in kmerfile:
	if line=='': break
	line=line.strip()
	if line.startswith("kmer"): continue
	split=line.split()
	kmer=split[0]
	chrom=split[1]
	pos=int(split[2])
	if chrom not in cdict: cdict[chrom]={}
	mybin=int(pos/1000)
	mybin=mybin*1000
	if mybin in cdict[chrom]:
		if kmer in cdict[chrom][mybin]: cdict[chrom][mybin][kmer]+=1
		if kmer not in cdict[chrom][mybin]: cdict[chrom][mybin][kmer]=1
	if mybin not in cdict[chrom]: cdict[chrom][mybin]={kmer:1}
kmerfile.close()
print "Counted all the kmers. Writing to new file."

mychroms=cdict.keys()
for c in sorted(mychroms):
	mywins=cdict[c].keys()
	for w in sorted(mywins):
		mykmers=cdict[c][w].keys()
		for k in sorted(mykmers):
			newfile.write(str(c)+'\t'+str(w)+'\t'+str(k)+'\t'+str(cdict[c][w][k])+'\n')

newfile.close()

