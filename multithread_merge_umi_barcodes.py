import sys
import multiprocessing as mp
from collections import Counter
from collections import defaultdict

sample_barcode_file = sys.argv[1] ## barcode file

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

def merge_umi(data):
	sample = data["sample"]
	filePath = sys.argv[2] + "/"
	inFile = open(filePath+sample+".aligned.txt",'r').read().splitlines()[::2]
	outFile = open(filePath+sample+".new.aligned.txt",'w')
	barcodeFile = open(filePath+sample+".newi.barcodes.txt",'w')

	bcds = []
	for line in inFile:
		fds = line.split("\t")
		n = fds[0]
		brcd = fds[1]
		seq = fds[2]
		mapped = fds[3]
		if mapped.count("X") < 10:
			outFile.write(brcd + "\t" + seq + "\t" + mapped + "\n" + "\n")
			bcds.append(brcd)
	bcounter = Counter(bcds)
	sbcds = sorted(bcounter.items(), key=lambda i: i[1], reverse=True)
	other = sorted(bcounter.items(), key=lambda i: i[1], reverse=True)
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
			barcodeFile.write(b[0] + "\t" + ";".join(sims[b[0]]) + "\t" + str(bcounter[b[0]]) + "\t" +str(b[1])+"\n")
	outFile.close()
	barcodeFile.close()

samples = []
jobs = []
with open(sample_barcode_file) as file:
	for line in file:
		line = line.split('\t')
		sample = line[0]
		samples.append(sample)
		job = {}
		job["sample"] = sample
		jobs.append(job)

pool = mp.Pool(processes=mp.cpu_count())
multithread = pool.map(merge_umi,jobs)