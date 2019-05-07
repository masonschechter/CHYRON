import pickle
import prefixTree as pTree
# from graphviz import Digraph, Source, nohtml, Graph
import itertools
import pandas as pd
import math
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
# from sklearn.cluster import hierarchical, KMeans
from scipy.spatial.distance import pdist, squareform
# from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
import operator
import numpy as np


# print(f"Number of concurrent experiments: {int(len(wells)/4)}")
# for exp in range(int(len(wells)/4)):
	# print(f"Wells in experiment {exp+1}: {wells[(exp*4):(exp*4)+4]}")

with open('darkgreenLT_wellDict.pkl', 'rb') as file:
	wellDict = pickle.load(file)


wells = []
for x in range(1,17):
	wells.append(f"well_{str(x)}")


pooledTree = pTree.PrefixTree(treeFile='darkgreenLT_pooled_insTree.pkl')
pooledNodes = [x for x in pooledTree.nodes]


binaryDF = pd.DataFrame(index=[pooledNodes])




for well in wells:
	wellTree = pTree.PrefixTree(treeFile=f"darkgreenLT_{well}_insTree.pkl")
	wellNodes = [x for x in wellTree.nodes]
	wellVector = []
	for n in pooledNodes:
		if n in wellNodes:
			wellVector.append(1)
		else:
			wellVector.append(0)
	binaryDF[well] = wellVector

binaryDF = binaryDF.transpose()
# print(binaryDF)
dist = pdist(binaryDF, 'jaccard')
# z = hierarchy.linkage(dist, 'average', optimal_ordering=True)
# plt.figure()

# dn = hierarchy.dendrogram(z)
# plt.show()
# z = hierarchy.fclusterdata(z, 8)
# plt.figure()

# dn = hierarchy.dendrogram(z)
# plt.show()


#### Generating insertion trees ####

conWells = []
for x, y in itertools.combinations(wells, 2):
	conWells.append(f"well_{x[5:]}_{y[5:]}")

# print(f"Number of concurrent experiments: {int(len(wells)/4)}")
# for exp in range(int(len(wells)/4)):
	# print(f"Wells in experiment {exp+1}: {wells[(exp*4):(exp*4)+4]}")

distDF = pd.DataFrame(columns=[x for x in wells], index=[x for x in wells])
# distDF = pd.DataFrame(columns=[x for x in conWells], index=[x for x in conWells])

wellsAndDistances = {}

for x, y in itertools.combinations(wells, 2): ## wells
	tree1 = pTree.PrefixTree(treeFile=f"darkgreenLT_well_{x[5:]}_insTree.pkl")
	tree2 = pTree.PrefixTree(treeFile=f"darkgreenLT_well_{y[5:]}_insTree.pkl")
	# print(f"Comparing {x} and {y}")
	consensusTree = pTree.PrefixTree()
	t1List = [x for x in tree1.nodes]
	t2List = [x for x in tree2.nodes]
	consensusNodes = list(set(t1List).intersection(t2List))
	consensusNodes.sort()
	for n in consensusNodes:
		# print(f"Inserting {n}")
		node = pTree.Node(n,None)
		consensusTree.insert(node,consensusTree.getRoot())
	print(f"Well {x[5:]} has {len(tree1.nodes)} nodes")
	print(f"Well {y[5:]} has {len(tree2.nodes)} nodes")
	print(f"Consensus tree for {x} and {y} has {len(consensusTree.nodes)} nodes")
	# consensusTree.export(f"well_{x[5:]}_{y[5:]}_insTree.pkl")
	con = len(consensusTree.nodes)
	lenX = len(tree1.nodes)
	lenY = len(tree2.nodes)
	# con = weight(consensusNodes)
	# list_x = weight(t1List)
	# list_y = weight(t2List)

# # 	##### Kulczynski-2 Similarity
	# K = True
	# J = False
	# O = False
	# distDF.at[x, y] = .5*((con/(con+lenX)) + (con/(con+lenY)))
	# distDF.at[y, x] = .5*((con/(con+lenX)) + (con/(con+lenY)))
	# distDF.at[x, x] = 1
	# distDF.at[y, y] = 1
# # 	##### Jaccard Similarity (J)
	K = False
	J = True
	O = False
	Jaccard = con/(lenX + lenY - con)
	distDF.at[x, y] = Jaccard
	distDF.at[y, x] = Jaccard
	distDF.at[x, x] = 1
	distDF.at[y, y] = 1
	wellsAndDistances[f'{x}-{y}'] = Jaccard
	orderedPairs = sorted(wellsAndDistances.items(), key=operator.itemgetter(1), reverse=True)


# 	####Ochiai Similarity
	# K = False
	# J = False
	# O = True
	# term_1 = con/(con+lenX)
	# term_2 = con/(con+lenY)
	# distDF.at[x, y] = math.sqrt(term_1*term_2)
	# distDF.at[y, x] = math.sqrt(term_1*term_2)
	# distDF.at[x, x] = 1
	# distDF.at[y, y] = 1

if J:
	print("Jaccard Similarity:")
elif K:
	print("Kulczynski-2 Similarity")
elif O:
	print("Ochiai Similarity")
print(distDF)
distDF.to_csv("instree_Aug1-4_C50_L5.csv")

seenWells = []
final = []
stderror = {}
for wells,dist in orderedPairs:
		well1,well2 = wells.split('-')
		# print(f"Working on {wells}")
		col = [x for x in distDF.loc[well1] if x < 1]
		if len(col) > 1:
			std = np.std(col)
			mean = np.mean(col)
			error = (dist-mean)/std
			stderror[wells] = error
		else:
			stderror[wells] = 0
		if well1 not in seenWells and well2 not in seenWells:
			final.append(wells) # Wells is the string of both wells
			seenWells.append(well1)
			seenWells.append(well2)
			# distDF.loc[well1] = np.nan
			# distDF[well2] = np.nan
			# distDF.loc[well2] = np.nan
			# distDF[well1] = np.nan
		# print(f"Current Value: {dist} | Mean: {mean} | Standar Div: {std}")
		# print(f"We are {error} standard deviations away for {well1}\n")
print("**********Lowest Level***************")
totalStd = 0
for winnerWells in final:
	print(f"{winnerWells} with a standard err of {stderror[winnerWells]}")
	totalStd += stderror[winnerWells]
print(f"Average Standard Err: {round(totalStd / len(final),2)}")
print(f"Final confidence of {round(norm.cdf(totalStd / len(final)) * 100,2)}%")


# for item in orderedPairs: ##debug
# 	print(item)


z = hierarchy.linkage(distDF, 'average')#, optimal_ordering=True)
plt.figure()

dn = hierarchy.dendrogram(z)
plt.show()

# distDF.to_csv("Jaccard_7-15_count10.csv")


# print(distDF.idxmax())
# print(distDF.idxmax(axis=1))