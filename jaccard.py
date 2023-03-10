import pickle
import pandas as pd
from itertools import combinations, product
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.stats import norm
import operator
import numpy as np

with open('darkgreenLT_insDict.pkl', 'rb') as file:
	insDict = pickle.load(file)

LENGTH_CUTOFF = (7,16) ## inclusive, exclusive
COUNT_CUTOFF = 50
RECOVERY_EFFICIENCY = .05

def jaccard_calc(l1,l2):
	a,b,c,d = (0,0,0,0)
	for i in range(len(l1)):
		if l1[i] == l2[i]:
			if l1[i] == 1:
				a+=1
				continue
			else:
				d+=1
				continue
		else:
			if l1[i] == 1:
				b+=1
				continue
			else:
				c+=1
				continue
	# print(a,b,c,d)
	return (a)/(a+b+c)

def insertionFinder(wellPairs, vectors):
	for pair in wellPairs:
		with open(f"{pair[0]}_{pair[1]}_insertions.txt",'w') as file:
			insertions = vectors["wells"][pair[0]][0].intersection(vectors["wells"][pair[1]][0])
			for insertion in insertions:
				file.write(insertion+'\n')

def nextLevelInsertionFinder(wellPairs2, vectors):
	for pair in wellPairs2:
		with open(f"{pair[0]}_{pair[1]}_{pair[2]}_{pair[3]}_insertions.txt",'w') as file:
			set1 = vectors["wells"][pair[0]][0].intersection(vectors["wells"][pair[1]][0])
			set2 = vectors["wells"][pair[2]][0].intersection(vectors["wells"][pair[3]][0])
			insertions = set1.intersection(set2)
			for insertion in insertions:
				file.write(insertion+'\n')


# c =0
# for ins in insDict:
# 	print(ins, insDict[ins])
# 	c +=1
# 	if c >=10:
# 		break
# exit()

wells = []
vectors = {}
vectors["pool"] = [set(),[]]
vectors["wells"] = {}

group_1 = ["0","1","2","3","5","6","10","14"]
group_2 = ["4","7","8","9","11","12","13","15"]

for x in range(1,17):
	wells.append(f"well_{str(x)}")
	vectors["wells"][f"well_{str(x)}"] = [set(),[]]

# for i in range(1,5):
# 	vectors["wells"][f"well_{str(i)}"] = [set(),[]]
# 	wells.append(f"well_{str(i)}")
	# print(vectors["wells"])

# for _ in group_2:
# 	vectors["wells"][f"well_{_}"] = [set(),[]]
# 	wells.append(f"well_{_}")

wellPairs = [('well_7','well_16'), ('well_1','well_10'), ('well_2','well_11'), ('well_4','well_9'),
			('well_15','well_3'), ('well_12','well_14'), ('well_13','well_6'), ('well_8','well_5')]

wellPairs_2 = [('well_7','well_16', 'well_4','well_9'), ('well_2','well_11', 'well_1','well_10'),
			('well_15','well_3', 'well_12','well_14'), ('well_13','well_6', 'well_5','well_8')]

for ins in insDict:
	if not ins == 'ROOT':
		if len(ins) in range(*LENGTH_CUTOFF):
			for well in insDict[ins]["counts"]:
				if well in wells:
					if insDict[ins]["counts"][well] >= COUNT_CUTOFF:
						if np.random.choice([0,1], p=[1-RECOVERY_EFFICIENCY, RECOVERY_EFFICIENCY]):
							vectors["pool"][0].add(ins)
							vectors["wells"][well][0].add(ins)
				else:
					print("couldnt find well")
# nextLevelInsertionFinder(wellPairs_2, vectors)

distDF = pd.DataFrame(columns=[x for x in wells], index=[x for x in wells], dtype='float64')


for ins in vectors["pool"][0]:
	for well in vectors["wells"]:
		if ins in vectors["wells"][well][0]:
			vectors["wells"][well][1].append(1)
		else:
			vectors["wells"][well][1].append(0)

for well in vectors["wells"]:
	print(f'{well} length {len(vectors["wells"][well][0])}')

# for well in vectors["wells"]:
	# print(f'{well} length {len(vectors["wells"][well][1])}')

wellsAndDistances = {}
for x,y in combinations(wells, 2):
	d = jaccard_calc(vectors["wells"][x][1],vectors["wells"][y][1])
	distDF.at[x, y] = d #round(d,1)
	distDF.at[y, x] = d #round(d,1)
	distDF.at[y, y] = 1
	distDF.at[x, x] = 1
	wellsAndDistances[x+'-'+y] = d

orderedPairs = sorted(wellsAndDistances.items(), key=operator.itemgetter(1), reverse=True)
seenWells = []
final = []
stderror = {}

for wells,dist in orderedPairs:
		well1,well2 = wells.split('-')
		if well1 not in seenWells and well2 not in seenWells:
			final.append(wells) # Wells is the string of both wells
			seenWells.append(well1)
			seenWells.append(well2)
			col = [x for x in distDF.loc[well1] if x > 0]
			row = [y for y in distDF.loc[well2] if y > 0]
			row_col = []
			row_col.extend(row+col)
			# row_col.remove(dist)
			if len(col) > 1 and len(row) > 1:
				std = np.std(row_col)
				mean = np.mean(row_col)
				error = (dist-mean)/std
				stderror[wells] = error
			else:
				stderror[wells] = 0

print(distDF)
print("**********Lowest Level***************")
totalStd = 0
for winnerWells in final:
	print(f"{winnerWells} with a standard err of {stderror[winnerWells]}")
	totalStd += stderror[winnerWells]
print(f"Average Standard Err: {round(totalStd / (len(final)-1),2)}")
print(f"Final confidence of {round(norm.cdf(totalStd / len(final)) * 100,2)}%")
# if SYSTEM["SHOW_DENDROGRAM"]:
# distDF = distDF.replace(np.nan,1)
	# print(z)
	# print(dn)
	# plt.savefig("static/results.jpg")
	# return (returnString,distDF)
distDF.to_csv(f'jaccard_darkgreenLT_C{str(COUNT_CUTOFF)}_L{str(LENGTH_CUTOFF[0])}-{str(LENGTH_CUTOFF[1])}_Recovery_{str(int(RECOVERY_EFFICIENCY*100))}.csv')
z = hierarchy.linkage(distDF, 'average') ##there are a few clustering choices here. for UPGMA, a standard algorithm, use 'average' instead of 'ward'
# plt.figure(figsize=(14,6),dpi=100)
dn = hierarchy.dendrogram(z, labels=distDF.index, leaf_rotation=90)
plt.show()
