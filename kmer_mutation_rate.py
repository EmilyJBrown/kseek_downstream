## take a table of normalized kmer counts, a table of groups, and a divergence times, calculate mutation rate in kmers

import sys

kcountsfh=sys.argv[1]
popsfh=sys.argv[2]
divtime=float(sys.argv[3])
newfh=sys.argv[4]

kcounts=open(kcountsfh, 'r')
pops=open(popsfh, 'r')
newfile=open(newfh, 'w')

popdict={}

for line in pops:
	if line=='': break
	line=line.strip()
	split=line.split()
	lineid=split[0]
	mypop=split[1]
	if lineid not in popdict: popdict[lineid]=mypop

linedict={}

for line in kcounts:
	B10count=0
	B10mean=0
	if line=='': break
	line=line.strip()
	split=line.split()
	if line.startswith("lines"):
		newfile.write(str("kmer")+'\t'+str("length")+'\t'+str("mean"))
		for i in range(1,len(split)):
			if i not in linedict: linedict[i]=split[i]
			if split[i] in popdict:
				newfile.write('\t'+str(split[i]))
		newfile.write('\n')
		continue
	kmer=split[0]
	mylen=kmer.split('/')[0]
	mylen=len(mylen)
	newfile.write(str(kmer)+'\t'+str(mylen))
	for i in range(1,len(split)):
		myline=linedict[i]
		if myline not in popdict: continue
		if myline in popdict:
			B10mean+=float(split[i])
			B10count+=1
	B10mean=float(B10mean)/float(B10count)
	newfile.write('\t'+str(B10mean))
	for i in range(1, len(split)):
		myline=linedict[i]
		if myline not in popdict: continue
		if myline in popdict:
			mycount=float(split[i])
			gen=float(divtime)
			mymu=float(mycount-B10mean)/float(gen)
			newfile.write('\t'+str(mymu))
	newfile.write('\n')

newfile.close()
kcounts.close()
pops.close()
