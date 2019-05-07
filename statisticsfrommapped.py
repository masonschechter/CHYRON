import sys
import const
import os
import re
import numpy as np
from matplotlib import pyplot as plt
import collections
from collections import defaultdict
import regex
import os.path

'''
sample_numbers
directory of aligned files
reference file
output file directory and base name
insertion_dir
'''


bads1 = const.bads1
bads2 = const.bads2
R_H = const.R_H
MAX_INSDEL = const.MAXINSDEL

## A File with the reference seq inside
ref = open(sys.argv[3],'r').read().splitlines()[0]





acgt = {"A":"T" , "C":"G", "G":"C", "T":"A"}
## Sample Numbers
sample_files = open(sys.argv[1],'r').read().splitlines()

#print sample_files
sample_nums = []
cut_dict = {}
if const.FIND_CUTSIDE:
	## Name \tab 21nt sequence around cut side
	## I don't try to find the cut side anymore; it looks
	## The forward sequence is the primary
	for line in sample_files:
		fds = line.split('\t')
		whatwewant = fds[3].upper()
		if whatwewant in ref:
			this_cut = ref.index(whatwewant)+9
			cut_dict[fds[0]] = this_cut
			sample_nums.append(fds[0])
		else:
			whatwewant = "".join([acgt[i] for i in fds[3].upper()[::-1]])
			if whatwewant in ref:
				this_cut = ref.index(whatwewant)+9
				cut_dict[fds[0]] = this_cut
				sample_nums.append(fds[0])
else:
	for line in sample_files:
		fds = line.split()
		sample_nums.append(fds[0])
		cut_dict[fds[0]] = const.CUT
	
#print cut_dict
## directory to *.aligned.txt
## make sure about the *.aligned.txt 
ALIGN_DIR = sys.argv[2]
INS_DIR = sys.argv[5]


## Pre Name
out_dir = sys.argv[4]
outfile = open(sys.argv[4] + "-stats.txt",'w')
nucfile = open(sys.argv[4] + "-nucs.txt" ,'w')
extnucfile = open(sys.argv[4] + "-extnucs.txt", 'w')
size_counter_file = open(sys.argv[4] + "-sizecounters.txt", 'w')
start_counter_file = open(sys.argv[4] + "-startcounters.txt" ,'w')


def prelength(s):
	count = 0
	for l in s:
		if l=="I":
			count+=1
	return count

def prelengthbyread(s):
	count = 0
	for l in s:
		if l=="D":
			count+=1
	return count



nuchdrs = "\t".join(["Sample_number", "Ins_A", "Ins_C", "Ins_G", "Ins_T"])
nucfile.write(nuchdrs + "\n")


#stathdrs = "\t".join(["AroundCutPureInsertion", "AroundCutPureDeletion","AroundCutMutations","AroundCutInsertions","AroundCutDeletion(Regions)","AroundCut-oneIns-oneMut-noDel", "AroundCut-oneIns-twoMut-noDel", "AroundCut-oneIns-moreMut-noDel"])
stathdrs = "\t".join(["AroundCutPureInsertion", "AroundCutPureDeletion","AroundCutMutations","NoInsNoDelAroundCut","NoInsNoDelNoMutAroundCut"])

outfile.write("SAMPLE" + "\t" + "NumOfReads" + "\t" + stathdrs + "\n")

for this_sample in sample_nums:
	print "stats from " + this_sample
	CUT = cut_dict[this_sample]

	if os.path.exists(ALIGN_DIR+"/"+this_sample+".aligned.txt"):
		infile = open(ALIGN_DIR+"/"+this_sample+".aligned.txt",'r').read().splitlines()
		## This is to keep track of insertions
		newfile = open(INS_DIR+"/"+this_sample+".newinfo.txt", 'w')
		newfileupdated = open(INS_DIR+"/"+this_sample+".newinfo.withlengths.txt", 'w')

		print " working on stats " + this_sample
		mapped_seqs = []
		init_seqs = []
		#####################################################################
		## Just reading the file and parsing sequence and aligned MXDI seq
		#####################################################################
		for line in infile[0::2]:
			seq = line.split("\t")[1]
			if seq.count("X") <= 20   and  seq.count("D") < 0.5*len(ref):
				init_seqs.append(line.split('\t')[0])
				mapped_seqs.append(seq)
		print str(len(mapped_seqs)) + " reads are in this sample"
		print str(len(init_seqs)) + " reads are in this sample"

		####################################################################
		## Initialization for each sample
		#####################################################################
		inocc = 0.0      ##Number of seqs with any insertion
		delocc= 0.0      ##Number of seqs with any deletion
		xocc = 0.0       ##Number of seqs with any mutation
		inanddel =0.0    ##Number of seqs with deletion and insertion
		morein = 0.0     ##Number of seqs with more than one region of insertion
		moredel = 0.0    ##Number of seqs with more than one region of deletion
		morex = 0.0      ##Number of seqs with more than one base of mutation
		xdel = 0.0       ##Number of seqs with mutation and deletion
		xi = 0.0         ##Number of seqs with mutation and insertion
		x2del = 0.0      ##Number of seqs with more than one mutation and deletion
		x2i = 0.0        ##Number of seqs with more than one mutation and insertion
		wt = 0
		wtx = 0
		pureinsertion_start = []
		pureinsertion_size = []
		puredeletion_start = []
		puredeletion_size = []
		pureinsertionnucs = []
		puredeletionnucs = []
		pureinsertionnucsdict = defaultdict(list)
		puredeletionnucsdict = defaultdict(list)
		pureinsertioncounter = defaultdict(list)
		insertion_start = []
		insertion_size = []
		deletion_start = []
		deletion_size = []
		insertionnucs = []
		deletionnucs = []
		mutation_inits_start = []
		mutation_start = []
		xmorei = 0
		insertion_types = []

		for seqind in range(len(mapped_seqs)):
			seq = mapped_seqs[seqind]
			initseq = init_seqs[seqind]
			#####################################################################
			## Just finding the location as observed in the output of alignments
			D = re.compile('DD*')
			Dpos_inits = [[m.start(),m.end()] for m in D.finditer(seq)]
			I = re.compile('II*')
			Ipos_inits = [[m.start(),m.end()] for m in I.finditer(seq)]
			X = re.compile('X')
			Xpos_inits = [[m.start(),m.end()] for m in X.finditer(seq)]
			#####################################################################

			Xpos_inits_fixed = []
			Ipos_inits_fixed = []
			Dpos_inits_fixed = []

			Xpos_fixed = []
			Ipos_fixed = []
			Dpos_fixed = []
			#####################################################################
			## Get the locations based on reference sequence
			for p in Xpos_inits:
				Xpos_inits_fixed.append(p[0] - prelength(seq[:p[0]]))
			for p in Dpos_inits:
				Dpos_inits_fixed.append([p[0] - prelength(seq[:p[0]]), p[0] - prelength(seq[:p[0]]) + p[1] - p[0]] )
			## This one is not actually right, it just shows the length not the right end on reference
			for p in Ipos_inits:
				Ipos_inits_fixed.append([p[0] - prelength(seq[:p[0]]), p[0] - prelength(seq[:p[0]]) + p[1] - p[0]] )
			#####################################################################
			## In case we have a cut side
			if True:
				for mut_occ in Xpos_inits_fixed:
					if mut_occ < (CUT+R_H+1) and mut_occ>= (CUT-R_H) :
						Xpos_fixed.append(mut_occ)

				for insertion_occ in Ipos_inits_fixed:
					if True:
					#if len(Dpos_inits_fixed) ==0 and len(Ipos_inits_fixed)==1 and len(Xpos_fixed) == 0 :
						if insertion_occ[0] < (CUT+R_H+1) and insertion_occ[0] >= (CUT-R_H) :
							if insertion_occ[1] - insertion_occ[0] <= MAX_INSDEL:
								Ipos_fixed.append(insertion_occ)
				 
				for deletion_occ in Dpos_inits_fixed:
					if True:
					#if len(Ipos_inits_fixed) ==0 and len(Dpos_inits_fixed)==1 and len(Xpos_inits_fixed) == 0 :
						if (deletion_occ[0] < (CUT+R_H+1) and deletion_occ[0] >= (CUT-R_H)) or \
								(deletion_occ[1] <= (CUT+R_H+1) and deletion_occ[1] > (CUT-R_H)) or \
								(deletion_occ[0] < (CUT-R_H) and deletion_occ[1] > (CUT+R_H+1)):
							Dpos_fixed.append(deletion_occ)
						elif (deletion_occ[0] < (CUT-R_H) and deletion_occ[1] > (CUT+R_H+1)) or \
						(deletion_occ[0] in range(CUT-R_H,CUT+R_H+1) and deletion_occ[1] > (CUT+R_H+1)) or \
						(deletion_occ[0] < (CUT-R_H) and deletion_occ[1] in range(CUT-R_H,CUT+R_H+2)):
							if deletion_occ[1] - deletion_occ[0] <= MAX_INSDEL:
								Dpos_fixed.append(deletion_occ)

			## get the location of Insertions in this seq
			Ipos_all_second = []
			Ipos_second = [] 
			for p in Ipos_inits:
				Ipos_all_second.append([p[0] - prelengthbyread(seq[:p[0]]), p[0] - prelengthbyread(seq[:p[0]]) + p[1] - p[0]])
				if [p[0] - prelength(seq[:p[0]]), p[0] - prelength(seq[:p[0]]) + p[1] - p[0]] in Ipos_fixed:
					Ipos_second.append([p[0] - prelengthbyread(seq[:p[0]]), p[0] -prelengthbyread(seq[:p[0]]) + p[1] - p[0]]) 



			##########################################################
			##To get simple counts of base pairs and regions
			## Im not using most of these parameters
			#########################################################
			## Anywhere
			muts_init = len(Xpos_inits_fixed)
			## Around cut
			muts = len(Xpos_fixed)

			## Anywhere
			mutation_inits_start+=Xpos_inits_fixed
			## Around cut
			mutation_start+=Xpos_fixed
			
			#############################################################
			# Mutation, yes or no
			if muts > 0:
				xocc+=1
			# One mutation any where, one mutation and no deletion around cut              
			if muts == 1 and len(Ipos_fixed) ==1 and len(Dpos_fixed) == 0:
				xi+=1            
			# Two mutations any where, one mutation and no deletion around cut  
			if muts == 2 and len(Ipos_fixed) ==1 and len(Dpos_fixed) == 0:
				x2i+=1
			# More than two mutations any where, one mutation and no deletion around cut  
			if muts > 2 and len(Ipos_fixed) ==1 and len(Dpos_fixed) == 0 :
				xmorei+=1
			# Only one Insertion around Cut, nothing else around Cut
			if len(Ipos_fixed) == 1 and len(Dpos_fixed) == 0 and len(Xpos_fixed) == 0 :
				inocc+=1
			# Insertion anywhere
			if len(Ipos_inits_fixed) > 0:
				morein+=1
			# Only one Deletion around Cut, nothing else around Cut        
			if len(Dpos_fixed) == 1 and len(Ipos_fixed) == 0 and len(Xpos_fixed) == 0 :
				delocc+=1
			# Deletion anywhere
			if len(Dpos_inits_fixed) > 0 :
				moredel +=1
			# Insertion and Deletion around Cut
			if len(Ipos_fixed) > 0 and len(Dpos_fixed)>0 :
				inanddel+=1
			# No Insertion No deletion around Cut
			if len(Ipos_fixed) == 0 and len(Dpos_fixed) == 0:
				wt+=1
				if len(Xpos_fixed) == 0:
					wtx+=1
					### This is to keep root
					newfile.write(initseq + "\t" + seq + "\t" + "-1" + "\t" + "ROOT" + "\n")
			#############################################################
			##For nucleotide probabilities, size and starting position
			#############################################################

			for i in Ipos_second:
				if len(Dpos_fixed) == 0 and muts == 0 and len(Ipos_fixed)==1 :
					#print i
					this_insertion = initseq[i[0]: i[1]]

		# 			#pureinsertionnucs+=this_insertion ## moved to each if block

					if len(this_insertion)>0 and len(this_insertion)<=12:
						insertion_types.append(this_insertion)#
						pureinsertionnucs+=this_insertion#
						pureinsertionnucsdict[min(i[1]-i[0],20)] += this_insertion
						pureinsertioncounter[min(i[1]-i[0],20)].append("this")
						pureinsertion_size.append(i[1] - i[0])
						pureinsertion_start.append(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]]))
						newfile.write(initseq + "\t" + seq + "\t" + str(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]])) + "\t" + initseq[i[0]: i[1]] + "\n")
						insertionnucs+=initseq[i[0]:i[1]]
						insertion_size.append(i[1] - i[0])
						insertion_start.append(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]]))
					elif (len(this_insertion)>12 and len(regex.findall("(%s){s<=2}"%this_insertion,bads1.upper()))==0 and len(regex.findall("(%s){s<=2}"%this_insertion, bads2.upper()))) == 0:
						if len(this_insertion)>=16:
							continue
							# with open('hg38.fa') as genome: ## check if insertion is a homologous recombination product from 209, 405, or human genome
							# 	for line in genome:
							# 		if not regex.findall("(%s){s<=3}"%this_insertion.upper(),line.upper()): ## no matches found in line
							# 			continue
							# 		else:
							# 			print "found HR products: " + this_insertion + ", in " + this_sample
							# 			break
								# else: ## no matches found in genome
								# 	insertion_types.append(this_insertion) #
								# 	pureinsertionnucs+=this_insertion
								# 	pureinsertionnucsdict[min(i[1]-i[0],20)] += this_insertion
								# 	pureinsertioncounter[min(i[1]-i[0],20)].append("this")
								# 	pureinsertion_size.append(i[1] - i[0])
								# 	pureinsertion_start.append(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]]))
								# 	newfile.write(initseq + "\t" + seq + "\t" + str(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]])) + "\t" + initseq[i[0]: i[1]] + "\n")
								# 	insertionnucs+=initseq[i[0]:i[1]]
								# 	insertion_size.append(i[1] - i[0])
								# 	insertion_start.append(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]]))
						else: # 12 < len(ins) < 16
							insertion_types.append(this_insertion) #
							pureinsertionnucs+=this_insertion
							pureinsertionnucsdict[min(i[1]-i[0],20)] += this_insertion
							pureinsertioncounter[min(i[1]-i[0],20)].append("this")
							pureinsertion_size.append(i[1] - i[0])
							pureinsertion_start.append(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]]))
							newfile.write(initseq + "\t" + seq + "\t" + str(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]])) + "\t" + initseq[i[0]: i[1]] + "\n")
							insertionnucs+=initseq[i[0]:i[1]]
							insertion_size.append(i[1] - i[0])
							insertion_start.append(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]]))
				# insertionnucs+=initseq[i[0]:i[1]]
				# insertion_size.append(i[1] - i[0])
				# insertion_start.append(i[0]+ prelengthbyread(seq[:i[0]])- prelength(seq[:i[0]]))
				

		#####################################################################
		###Looking at pure insertion
		insnum = len(pureinsertionnucs)/100.0
		if insnum == 0:
			insnum = 1
		purenuctemp = [pureinsertionnucs.count("A")/insnum, pureinsertionnucs.count("C")/insnum, pureinsertionnucs.count("G")/insnum, pureinsertionnucs.count("T")/insnum]

		####################################################################
		###Looking at any insertion 
		insnum = len(insertionnucs)/100.0
		if insnum == 0:
			insnum = 1

		nuctemp = [insertionnucs.count("A")/insnum, insertionnucs.count("C")/insnum, insertionnucs.count("G")/insnum, insertionnucs.count("T")/insnum]

		nucfile.write(this_sample + "\t"  + "\t".join(map(str,purenuctemp))  + "\t" + "====" + "\t" + "\t".join(map(str,nuctemp))  + "\n")

		##################################################################
		normed = len(mapped_seqs)/100.0
		if normed == 0 :
			normed = 1.0
		temp_insertion_types_counter = collections.Counter(insertion_types)
		insertion_types_counter = {}
		for m in temp_insertion_types_counter:
			if temp_insertion_types_counter[m] > 0:
				insertion_types_counter[m] = temp_insertion_types_counter[m]*100.0/len(pureinsertion_size)
				
		
		stattemp = "\t".join(map(str,[round(inocc/normed,2), round(delocc/normed,2),round(xocc/normed,2), round(wt/normed,2), round(wtx/normed,2), str(sum([len(i) for i in insertion_types])/float(max(1,len(insertion_types)))),str(len(temp_insertion_types_counter)), str(temp_insertion_types_counter)]))

		####################################################################
		##To keep Some details on distribution
		####################################################################
		pureinsertion_start_counter = collections.Counter(pureinsertion_start)
		pureinsertion_size_counter = collections.Counter(pureinsertion_size)

		insertion_inits_start_counter = collections.Counter(insertion_start)
		insertion_inits_size_counter = collections.Counter(insertion_size)

		'''
		plt.hist(insertion_inits_start_counter.keys(), weights=insertion_inits_start_counter.values(), bins=range(0,250), facecolor='grey', alpha=0.55)
		plt.hist(pureinsertion_start_counter.keys(), weights=pureinsertion_start_counter.values(), bins=range(0,250), facecolor='teal', alpha=0.75)
		plt.ylabel("Frequency")
		plt.title(this_sample+"_start")
		plt.savefig(sys.argv[4] + "-" + str(this_sample) + "-" + "start.png")

		plt.gcf().clear()
		plt.hist(insertion_inits_size_counter.keys(), weights=insertion_inits_size_counter.values(), bins=range(1,30), facecolor='grey', alpha=0.55)
		plt.hist(pureinsertion_size_counter.keys(), weights=pureinsertion_size_counter.values(), bins=range(1,30), facecolor='red', alpha=0.75)
		plt.ylabel("Frequency")
		plt.title(this_sample+"_size")
		plt.savefig(sys.argv[4] + "-" + str(this_sample) + "-" + "size.png")
		plt.gcf().clear()
		'''

		size_counter_file.write(this_sample+ "\t" + str(pureinsertion_size_counter) +   "\n")
		start_counter_file.write(this_sample+ "\t" + str(pureinsertion_start_counter) +   "\n")

		####################################################################

		outfile.write(this_sample + "\t" + str(len(mapped_seqs)) + "\t" + stattemp + "\n")
		print this_sample + "  is done"
		for this_size in range(1,21):
			insnum = len(pureinsertionnucsdict[this_size])/100.0
			if insnum == 0:
				insnum = 1
			extnuctemp = [pureinsertionnucsdict[this_size].count("A")/insnum, pureinsertionnucsdict[this_size].count("C")/insnum, pureinsertionnucsdict[this_size].count("G")/insnum, pureinsertionnucsdict[this_size].count("T")/insnum]
			this_line = this_sample + "\t" + str(this_size) + "\t" + str(len(pureinsertionnucsdict[this_size])) +"\t" + str(len(pureinsertioncounter[this_size])) + "\t" + "\t".join(map(str,extnuctemp)) +"\n"
			extnucfile.write(this_line)
		if len(pureinsertion_size) > 0 :
			extnucfile.write("Average_Insertion_Length" +"\t" + str(float(sum(pureinsertion_size))/len(pureinsertion_size)) + "\n")
		else:
			 extnucfile.write("Average_Insertion_Length\tNA" +"\n")
		extnucfile.write("\n")
		newfile.close()
		
		newfile_reopened = open(INS_DIR + "/"+this_sample+".newinfo.txt", 'r').read().splitlines()
		keep_track= []
		keep_length = {}
		for line in newfile_reopened:
			fds = line.split()
			lengg = len(fds[-1])
			iddd = fds[-2] + "=" + fds[-1]
			keep_length[iddd] = lengg
			keep_track.append(iddd)
		
		
		counts = collections.Counter(keep_track)
		nrmz = len(newfile_reopened)
		for i in counts:
			if keep_length[i] < 40:
				newfileupdated.write(i.split("=")[0] + "\t" + i.split("=")[1] + "\t" + str(keep_length[i]) + "\t" +str(counts[i]) + "\t" + str(counts[i]*100.0/nrmz) + "\n")

		newfileupdated.close()
		newcommand = "sort -k3,3rn -k4,4rn " + INS_DIR + "/" +this_sample + ".newinfo.withlengths.txt > "+ INS_DIR + "/"  + this_sample + ".newinfo.withlengths.sorted.txt"
		os.system(newcommand)
		
outfile.close()
nucfile.close()
size_counter_file.close()
start_counter_file.close()
extnucfile.close()
