import pickle
import sys
import pprint
import re

# with open('green_LT_insDict.pkl', 'rb') as pkl:
# 	insDict = pickle.load(pkl)

# insDict = {
	# insertion:{		
# 		"counts": {
#			"well" : int()
#		}
# 	}
# }
## sysargv[1] = "exp"_well_data_file.txt

insDict = {}
with open(sys.argv[1]) as samples:
	for line in samples:
		counter = 0
		c = 0
		count = 0
		insCount = 0
		line = line.strip().split(' ')
		withlengths = line[0]
		well_regex = re.search('(?<=_).*?(?=[.])', line[0]) ## finds well number between '_' and '.' of file name
		well = f"well_{str(int(well_regex.group())+1)}" 
		print(f"Working on {well}")

		with open(withlengths) as lengths:
			for line in lengths:
				insInfo = line.strip().split('\t')
				insertion = insInfo[1]
				insLen = insInfo[2]
				insCount = int(insInfo[3])
				insPerc = insInfo[4]
				if insertion == 'ROOT':
					print("found one")
					continue
				if insertion not in insDict: ## have we seen this insertion?
					insDict[insertion] = {"counts":{well:insCount}}
					c += 1
					continue
				if well in insDict[insertion]["counts"]: ## have we seen this insertion in this well?
					insDict[insertion]["counts"][well] += insCount
				else:
					insDict[insertion]["counts"][well] = insCount
			# print(f"{c} unique insertions in {well}")
			# print(f"{insCount} non-unique in {well}")

countDict = {}

for key in insDict:
	for w in insDict[key]["counts"]:
		if not w in countDict:
			countDict[w] = 1
		else:
			countDict[w] += 1


pp = pprint.PrettyPrinter(indent=4)
pp.pprint(countDict)

with open('darkgreenLT_insDict.pkl', 'wb') as file:
	pickle.dump(insDict, file)
