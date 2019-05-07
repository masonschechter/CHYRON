import sys
import subprocess
import glob
import multiprocessing as mp

refseqdir = open(sys.argv[1] + "_ref.txt", 'r').read().splitlines()[0]
sample_barcode_file = sys.argv[2] ## barcode file
IN_DIR = sys.argv[3]
ALIGN_DIR = sys.argv[4]

acgt = {"A":"T", "C":"G", "G":"C", "T":"A"}


def DXMI(data):
	sample = data["sample"]
	fwd_barcode = data["fwd_barcode"]
	rev_barcode = data["rev_barcode"]
	# infile = open(IN_DIR + "/" + sample + "_forward_" + fwd_barcode + "_" + rev_barcode + ".fastq",'r').read().splitlines()
	infile = open(IN_DIR + "/" + sample + ".fastq.assembled.fastq",'r').read().splitlines()
	alignmentfile = open(ALIGN_DIR + "/" + sample + ".aligned.txt",'w')

	print " DXMI on sample " + sample
	mapped_seqs = []
	checked_seqs = {}
	for indd in range(len(infile[0::4])):
		this_seq = infile[4*indd+1].upper()
		if this_seq not in checked_seqs:
			proc = subprocess.Popen('./mapp ' + this_seq + '  ' + refseqdir, stdout=subprocess.PIPE, shell=True)
			this_aligned = proc.stdout.read()
			checked_seqs[this_seq] = this_aligned
		else:
			this_aligned = checked_seqs[this_seq]
		alignmentfile.write(this_seq + "\t" + this_aligned + "\n")
		mapped_seqs.append(this_aligned)
	alignmentfile.close()
	print " finished DXMI on " + sample

sample_barcode_pairs = []
with open(sample_barcode_file) as file:
	for line in file:
		line = line.split('\t')
		sample = line[0]
		fwd_barcode = line[1]
		rev_barcode = line[2]
		sample_barcode_pairs.append([sample, fwd_barcode, rev_barcode])

jobs = [] # create list of jobs to be multithreaded
for sample_barcode_pair in sample_barcode_pairs:
	job = {}
	job["sample"] = sample_barcode_pair[0]
	job["fwd_barcode"] = sample_barcode_pair[1]
	job["rev_barcode"] = sample_barcode_pair[2]
	jobs.append(job)

pool = mp.Pool(processes=mp.cpu_count())
multithread = pool.map(DXMI,jobs)