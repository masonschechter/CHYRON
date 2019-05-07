import pickle
import prefixTree as pTree
# from graphviz import Digraph, Source, nohtml, Graph
from operator import itemgetter

# g = Graph('g', filename='insTree.gv',node_attr={'shape': 'record', 'height': '.1','width':'1.0'},graph_attr={'ranksep':".5"})
# tree = pTree.PrefixTree() #empty tree

with open('darkgreenLT_insDict.pkl', 'rb') as file:
	insDict = pickle.load(file)

# insDict = {
	# insertion:{
# 		"seqs" : {
#			"well" : []
#		}
# 		"umis" : {
#			"well" : []
#		}
# 		"counts": {
#			"well" : int()
#		}
# 	}
# }

# wellDict = {
# 	"well": {
# 		"ins": {
# 			"umi_count": int()
# 			"insCount": int()
# 		}
# 	}
# }
wellDict = {}
for ins in insDict:
	if not ins == 'ROOT':
		if len(ins) in range(7,16):
			for well in insDict[ins]["counts"]:
				if insDict[ins]["counts"][well] >= 50:
					if well not in wellDict:
						wellDict[well] = {}
						wellDict[well][ins] = {"insCount":0}
					else:
						wellDict[well][ins] = {"insCount":0}
					wellDict[well][ins]["insCount"] = insDict[ins]["counts"][well]

with open('darkgreenLT_wellDict.pkl', 'wb') as file:
	pickle.dump(wellDict, file)


# group_1 = ["0","1","2","3","5","6","10","14"]
# group_2 = ["4","7","8","9","11","12","13","15"]

# for _ in group_2:
# 	wellDict[f"well_{_}"] = {}

# for ins in insDict:
# 	if not ins == 'ROOT':
# 		if len(ins) in range(5,16):
# 			for well in insDict[ins]["counts"]:
# 				if well in wellDict:
# 					if insDict[ins]["counts"][well] >= 50:
# 						wellDict[well][ins] = {"insCount":insDict[ins]["counts"][well]}

# with open('darkgreenLT_wellDict.pkl', 'rb') as file:
# 	wellDict = pickle.load(file)

# with open('darkgreenLT_wellDict.pkl', 'rb') as file:
# 	wellDict = pickle.load(file)
	
for well in wellDict:
	# print(f"{len(wellDict[well])} insertions for {well}, len 8-16")
	q = []
	# g = Graph('g', filename=f'{well}_insTree.gv',node_attr={'shape': 'record', 'height': '.1','width':'1.0'},graph_attr={'ranksep':".5"})
	tree = pTree.PrefixTree()
	for ins in wellDict[well]:
		if len(ins) in range(5,16):
			if wellDict[well][ins]["insCount"] >= 50:
				q.append(ins)
	q = sorted(q, key=itemgetter(0))
	print(f"{well} has {len(q)} insertions")
	for ins in q:
		# print(f"Inserting {ins}")
		node = pTree.Node(ins,None)
		tree.insert(node,tree.getRoot())
		# print(tree.nodes[node.id].parent)
	# for node in tree.nodes:
	# 	if tree.nodes[node].parent:
	# 		if tree.nodes[node].isGhost():
	# 			g.node(tree.nodes[node].id,_attributes={"shape":"point"})
	# 		else:
	# 			g.node(tree.nodes[node].id)
	# 		g.edge(tree.nodes[node].parent.id,tree.nodes[node].id,)#_attributes={'penwidth':str(tree.edges[tree.nodes[node].parent.id+"-"+tree.nodes[node].id]/2)})
	# 	else:
	# 		g.node(tree.nodes[node].id,_attributes={"shape":"point"})
	tree.export(f"darkgreenLT_{well}_insTree.pkl")
# counter = 0
# q = []
# for key in insDict:
# 	if counter >= 10000:
# 		break;
# 	if len(insDict[key]["ins"]) < 15:
# 		if insDict[key]["insCount"] > 3:
# 			q.append(insDict[key]["ins"])
# 	counter +=1
# q.sort();
# for key in q:
# 	print("Inserting: "+key)
# 	node = pTree.Node(key,None)
# 	tree.insert(node,tree.getRoot())
# 	counter+=1

wells = [x for x in wellDict]
pool = []
pooledTree = pTree.PrefixTree()
for well in wells:
	wellTree = pTree.PrefixTree(treeFile=f"darkgreenLT_{well}_insTree.pkl")
	wellNodes = [x for x in wellTree.nodes]
	for n in wellNodes:
		if n not in pool:
			pool.append(n)
for n in pool:
	node = pTree.Node(n,None)
	pooledTree.insert(node, pooledTree.getRoot())
pooledTree.export('darkgreenLT_pooled_insTree.pkl')
