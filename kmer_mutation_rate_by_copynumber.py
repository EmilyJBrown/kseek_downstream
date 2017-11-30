import sys

mutfh=sys.argv[1]
newfh=sys.argv[2]

mutfile=open(mutfh, 'r')
newfile=open(newfh, 'w')

for line in mutfile:
	if line=='': break
	line=line.strip()
	if line.startswith("kmer"):
		newfile.write(str(line)+'\n')
		continue
	split=line.split()
	mymean=float(split[2])
	newfile.write(str(split[0])+'\t'+str(split[1])+'\t'+str(split[2]))
	for i in range(3,len(split)):
		mymu=float(split[i])
		mynorm=float(mymu)/float(mymean)
		newfile.write('\t'+str(mynorm))
	newfile.write('\n')

mutfile.close()
newfile.close()