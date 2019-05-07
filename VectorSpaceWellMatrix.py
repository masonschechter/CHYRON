import numpy as np
from lib import loadFromWorkbook, loadFromPickle,distance_with_len,computeWellDistance
import pandas as pd
from itertools import combinations, product
from math import inf
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import sys, time
import multiprocessing as mp
import operator
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm

###############################
###### SYSTEM PARAMETERS ######
###############################
SYSTEM = {
"DECAY_RATE":2,
"MAGNITUDE": 1000,
"PREFIX_MATCH":1,
"LEN_CUTOFF":(9,16),
"COUNT_CUTOFF": 500,
"LOAD_FROM_EXCEL":False,
"EXCEL_DATA_FILE":"test for algorithm.xlsx",
# "PICKLE_DATA_FILE":"8_well_ins_dict.pkl",
"PICKLE_DATA_FILE":"darkgreenLT_insDict.pkl",
"SHOW_DENDROGRAM":False,
"MAX_THREADS":8
}
###############################
###############################
###############################
def main():
	print(f"**************************\nParameters:")
	for param in SYSTEM:
		print(f"{param}: {SYSTEM[param]}")
	print("**************************")
	if SYSTEM["LOAD_FROM_EXCEL"]:
		wellDict = loadFromWorkbook(SYSTEM["EXCEL_DATA_FILE"])
	else:
		wellDict = loadFromPickle(SYSTEM["PICKLE_DATA_FILE"],SYSTEM["LEN_CUTOFF"],SYSTEM["COUNT_CUTOFF"])

	wells = []
	for i in range(1,len(wellDict)+1):
		wells.append(f"well_{i}")
	for well in wells:
		print(f"{well} has {len(wellDict[well])} insertions")
	# exit()

	distDF = pd.DataFrame(columns=[x for x in wells], index=[x for x in wells])
	distDF = distDF.fillna(0.0)
	jobs = [] # Create a list of jobs that will be threaded.
	for wellX, wellY in combinations(wells, 2):
		# if wellX == wellY:
		# 	continue
		job = {}
		job["wellX"]=list(wellDict[wellX].keys())
		job["wellY"]=list(wellDict[wellY].keys())
		job["nameX"]=wellX
		job["nameY"]=wellY
		jobs.append(job)
		job2 = {}
		job2["wellX"]=list(wellDict[wellY].keys())
		job2["wellY"]=list(wellDict[wellX].keys())
		job2["nameX"]=wellY
		job2["nameY"]=wellX
		jobs.append(job2)

	#Initalize the pool
	if SYSTEM["MAX_THREADS"] > mp.cpu_count():
		print(f"Cannot satisfy config requirement of {SYSTEM['MAX_THREADS']} cores. This system only has {mp.cpu_count()} cores available")
		print(f"Using system max of {mp.cpu_count()} cores instead\n")
		SYSTEM["MAX_THREADS"] = mp.cpu_count()
	print(f"Running {len(jobs)} iterations")
	# result = computeWellDistance(jobs[0]) # for single result (debugging)
	# exit()
	pool = mp.Pool(processes=SYSTEM["MAX_THREADS"])
	multiResult = pool.map(computeWellDistance,jobs) #run it.
	#multiResult is a list of dictionaries all returned from the jobs
	wellsAndDistances = {} #Key is well pairing ('well4-well5'), and value is distance (score)
	for r in multiResult:
		x = r["nameX"]
		y = r["nameY"]
		result = r["distance"]
		name = x+'-'+y
		distDF.at[x, y] = result
		# distDF.at[y, x] = result
		distDF.at[x, x] = np.nan
		distDF.at[y, y] = np.nan
		if distDF.at[y, x] > 0:
			averageD = (distDF.at[y, x]+distDF.at[x,y]) / 2
			distDF.at[y, x] = averageD
			distDF.at[x, y] = averageD
			wellsAndDistances[y+'-'+x] = averageD
			wellsAndDistances[name] = averageD
	#Sort the combinations of well pairs by their distances
	orderedPairs = sorted(wellsAndDistances.items(), key=operator.itemgetter(1), reverse=True)
	for item in orderedPairs: # debug
		print(item) # debugging
	seenWells = []  # Keep track of what wells we've chosen so we don't choose it again
	final = []		# Final list of well pairs
	stderror = {}
	print(distDF)
	# distDF.to_csv('avgPrefix_Aug1-4_C50_L5.csv')
	#Now we go through the list IN ORDER, and take the first well combo we haven't seen yet (all the combinations)
	for wells,dist in orderedPairs:
		well1,well2 = wells.split('-')
		# print(f"Working on {wells}")
		col = [x for x in distDF.loc[well1] if x > 0]
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
	print(f"Average Standard Err: {round(totalStd / (len(final)-1),2)}")
	print(f"Final confidence of {round(norm.cdf(totalStd / len(final)) * 100,2)}%")
	
	if SYSTEM["SHOW_DENDROGRAM"]:
		distDF = distDF.replace(0,1)
		z = hierarchy.linkage(distDF, 'average') ##there are a few clustering choices here. for UPGMA, a standard algorithm, use 'average' instead of 'ward'
		plt.figure()
		dn = hierarchy.dendrogram(z)
		plt.show()


if __name__ == '__main__':
	startTime = time.time()
	main()
	print(f"Final Execution time: {time.time() - startTime} seconds. System had {mp.cpu_count()} cores available")