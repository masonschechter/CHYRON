import sys



def compare(seq1,seq2):
	counter =0
	for i in range(len(seq1)):
		if seq1[i]!=seq2[i]:
			counter+=1
		if counter>=2:
			return False
	if counter <=1:
		return True
	return False

name  = sys.argv[1]
infile = open(sys.argv[2]+"/"+name+".aligned.txt",'r').read().splitlines()[::2]

outfile = open(sys.argv[3]+"/"+name+".new.aligned.txt",'w')
barcodefile = open(sys.argv[3]+"/"+name+".newi.barcodes.txt",'w')

bcds = []
for line in infile:
	fds = line.split("\t")
	n = fds[0]
	brcd = fds[1]
	seq = fds[2]
	mapped = fds[3]
	if mapped.count("X") < 10:
		outfile.write(brcd + "\t" + seq + "\t" + mapped + "\n" + "\n")
		bcds.append(brcd)


from collections import Counter
bcounter = Counter(bcds)
sbcds = sorted(bcounter.items(), key=lambda i: i[1], reverse=True)
other = sorted(bcounter.items(), key=lambda i: i[1], reverse=True)
from collections import defaultdict
sims = defaultdict(list)
blacks = set()
for i in range(len(sbcds)-1):
	if i%500 == 0:
		print i
	if sbcds[i][0] not in blacks:
		this_count = sbcds[i][1]
		for j in range(i+1, len(sbcds)):
			if sbcds[j][0] not in blacks:
				if compare(sbcds[i][0],sbcds[j][0]):
					sims[sbcds[i][0]].append(sbcds[j][0])
					this_count+=sbcds[j][1]
					blacks.add(sbcds[j][0])
		sbcds[i]= (sbcds[i][0], this_count)

for b in sbcds:
	if b[0] not in blacks:
		barcodefile.write(b[0] + "\t" + ";".join(sims[b[0]]) + "\t" + str(bcounter[b[0]]) + "\t" +str(b[1])+"\n")

outfile.close()
barcodefile.close()

