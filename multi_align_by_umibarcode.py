import sys
from collections import Counter
from collections import defaultdict
import motility
import subprocess
import os

def make_pwm(list_of_seq):
	matrix = []
	#print list_of_seq[0]
	new_vertical = ["" for i in range(len(list_of_seq[0]))]
	for f in list_of_seq:
		for l in range(len(f)):
			new_vertical[l]+=f[l]


	for line in new_vertical:
		nts = line.upper().rstrip().replace("-", "N")
		freq = dict.fromkeys('ACGTN', 0.0)
		for nt in nts: freq[nt] += 1
		temp_list = [freq[bp] for  bp in "ACGT"]
		if max(temp_list) > 0.5*len(nts):	
	        	matrix.append(tuple((freq[bp] / len(nts)) for bp in "ACGT"))
	return motility.PWM(matrix)


def read_fasta(fp):
	name, seq = None, []
   	for line in fp:
   		line = line.rstrip()
        	if line.startswith(">"):
            		if name: yield (name, ''.join(seq))
            		name, seq = line, []
        	else:
            		seq.append(line)
    	if name: yield (name, ''.join(seq))

		    
infile = open(sys.argv[1],'r').read().splitlines()
similarfile = open(sys.argv[2],'r').read().splitlines()
similars = dict()
for line in similarfile:
	fds = line.split("\t")
	if len(fds) == 4:
		for other in fds[1].split(";"):
			similars[other] = fds[0]
	similars[fds[0]] = fds[0]

outfile = open(sys.argv[3] ,'w')
final_seqs = open(sys.argv[4],'w')

brcodes = []
seq_dict = defaultdict(list)
for line in infile[0::2]:
	fds = line.split("\t")
	this_one = fds[0]
	if this_one in similars:
		brcodes.append(similars[this_one])
		seq_dict[similars[this_one]].append(fds[1])

exp = sys.argv[5]
sample = sys.argv[6]
outdir = sys.argv[7]

import os
run = os.system 
proc = subprocess.Popen('module load mafft/7.305', shell=True)
c = 0
for bc in seq_dict:
	if c>-1 :
		if len(seq_dict[bc]) > 1 :
			this_list = seq_dict[bc]
			#print this_list
			in_temp = open(outdir + 'temp_for_umi/' + sample + '_temp_in.txt','w')
			for s in range(len(this_list)):
				in_temp.write(">" + str(s) + "\n" + this_list[s] + "\n")
			in_temp.close()
			out_temp = []

			#proc = subprocess.Popen('mafft --thread 16 --quiet --auto temp_in.txt > temp_out.txt', shell=True)
			#proc.communicate()
			#proc.wait()
			#proc.wait()
			result = run('mafft --thread 8 --quiet --auto ' + outdir + 'temp_for_umi/' + sample + '_temp_in.txt > '+ outdir + 'temp_for_umi/' + sample + '_temp_out.txt')
			with open(outdir + 'temp_for_umi/' + sample + '_temp_out.txt') as fp:
				for name, seq in read_fasta(fp):
					out_temp.append(seq.replace("-","n").upper())
			#os.remove("temp_in.txt")
			#os.remove("temp_out.txt")
			#if len(out_temp) == 0:
			#	print this_list
			#print out_temp[:5]
			pwm = make_pwm(out_temp)
			#pwm = motility.make_pwm(out_temp)
			max_score = pwm.max_score()
			best_sites = pwm.generate_sites_over(max_score)
			# ('GGAAACCGAA', 'GGAAACCGCA', 'GGAAACCGGA', 'GGAAACCGTA')
			#print best_sites
			#print motility.make_iupac_motif(best_sites)
			final_seqs.write(bc + "\t" + best_sites[0] + "\n")
			# GGAAACCGNA
			c+=1
			if c%25 == 0:
				print c
		else:
			c+=1
			final_seqs.write(bc + "\t" + seq_dict[bc][0] + "\n")
counts = Counter(brcodes)
counts_sorted = counts.most_common()
for f in counts_sorted:
	outfile.write(f[0] + "\t" + str(f[1]) + "\n")


outfile.close()
final_seqs.close()
