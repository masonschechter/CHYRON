import sys
import subprocess
import glob
'''
Input 1:
List of Samples

Input 2:
Directory to input merged no barcode reads

Input 3:
Directory to output alignments
'''
acgt = {"A":"T", "C":"G", "G":"C", "T":"A"}
##Numbers each in a separate line
sample_numbers=[sys.argv[2]]

##"$samplenumber"-merged.assembled.fastq are in this directory
IN_DIR = sys.argv[3]

##"$samplenumber".aligned.fastq go inside this directory
ALIGN_DIR = sys.argv[4]

FWD_BARCODE = [sys.argv[5]]
REV_BARCODE = [sys.argv[6]]

# refseqdir = open(sys.argv[1]+ '/' + sys.argv[1].split("/")[1]+ "_ref.txt", 'r').read().splitlines()[0]
refseqdir = open(sys.argv[1] + "_ref.txt", 'r').read().splitlines()[0]
#print refseqdir


for this_sample in sample_numbers:
	# infile = open(IN_DIR + "/" + this_sample + ".fastq.assembled.fastq",'r').read().splitlines()
	infile = open(IN_DIR + "/" + this_sample + "_forward_" + FWD_BARCODE[0] + "_" + REV_BARCODE[0] + ".fastq",'r').read().splitlines()
	#filename = glob.glob(IN_DIR+this_sample + '_reverse'+ '*.fastq')[0]
	#infile = open(filename,'r').read().splitlines()
	#subprocess.check_output(['ls','-l']) #all that is technically needed...
	alignmentfile = open(ALIGN_DIR+"/"+this_sample+".aligned.txt",'w')
	print " working on  sample " + this_sample
	mapped_seqs = []
	checked_seqs = {}
	for indd in range(len(infile[0::4])):
		this_seq = infile[4*indd+1].upper()
		#this_seq = "".join([acgt[i] for i in this_seq[::-1]])
		if this_seq not in checked_seqs:
			proc = subprocess.Popen('./mapp ' + this_seq + '  ' + refseqdir, stdout=subprocess.PIPE, shell=True)
			this_aligned = proc.stdout.read()
			checked_seqs[this_seq] = this_aligned
		else:
			this_aligned = checked_seqs[this_seq]
		alignmentfile.write(this_seq + "\t"+this_aligned + "\n")
		mapped_seqs.append(this_aligned)
	alignmentfile.close()

