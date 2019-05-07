import sys
import subprocess
'''
Input 1:
List of Samples

Input 2:
Directory to input merged no barcode reads

Input 3:
Directory to output alignments
'''

##Numbers each in a separate line
sample_numbers=[sys.argv[2]]

##"$samplenumber"-merged.assembled.fastq are in this directory
IN_DIR = sys.argv[3]

##"$samplenumber".aligned.fastq go inside this directory
ALIGN_DIR = sys.argv[4]

#### put this inside sample number loop vvvv ####
refseqdir = open(sys.argv[1] + "_ref.txt", 'r').read().splitlines()[0]
#refseqdir = open(sys.argv[1] + sample_numbers[0] + "_ref.txt", 'r').read().splitlines()[0]

for this_sample in sample_numbers:
	infile = open(IN_DIR+"/"+this_sample+".fastq.assembled.fastq",'r').read().splitlines()
	# refseqdir = open(sys.argv[1] + this_sample + "_ref.txt", 'r').read().splitlines()[0]
	#subprocess.check_output(['ls','-l']) #all that is technically needed...
	alignmentfile = open(ALIGN_DIR+"/"+this_sample+".aligned.txt",'w')
	print " map_togetDXMI on  sample " + this_sample
	mapped_seqs = []
	checked_seqs = {}
	for indd in range(len(infile[0::4])):
		this_one = infile[4*indd].split("_")[1]
		this_seq = infile[4*indd+1]
		if this_seq not in checked_seqs:
			proc = subprocess.Popen('./mapp ' + this_seq + " " + refseqdir, stdout=subprocess.PIPE, shell=True)
			this_aligned = proc.stdout.read()
			checked_seqs[this_seq] = this_aligned
		else:
			this_aligned = checked_seqs[this_seq]
		alignmentfile.write(infile[4*indd] +"\t" + this_one +  "\t" + this_seq + "\t"+this_aligned + "\n")
		mapped_seqs.append(this_aligned)
	alignmentfile.close()