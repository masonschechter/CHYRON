import sys

count_file = open(sys.argv[1],'r').read().splitlines()
mapped_file = open(sys.argv[2],'r').read().splitlines()


outfile = open(sys.argv[3],'w')
mapped_dict = {}
for line in mapped_file[::2]:
	fds = line.split("\t")
	mapped_dict[fds[-2]] = fds[-1]

notthere = 0.0
for line in count_file:
	fds = line.split('\t')
	m = fds[1]
	if m in mapped_dict:
		temp = m + "\t" + mapped_dict[m] + "\n"
		outfile.write(temp+"\n")
	else:
		notthere+=1
# print int(notthere/len(count_file))
outfile.close()
