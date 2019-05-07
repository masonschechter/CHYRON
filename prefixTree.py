import pickle

class Node():
	def __init__(self, name, parent, ghost=False):
		self.id = name
		self.parent = parent
		self.children = [];
		self.ghost = ghost

	def isGhost(self):
		return self.ghost

	def __str__(self):
		return self.id


class PrefixTree():
	def __init__(self,treeFile=None):
		if treeFile:
			self.loadFromFile(treeFile)
		else:
			# Create root node and add to dicts
			self.root = Node("0",None)
			self.edges = {"0":0}
			self.nodes = {"0":self.root}

	def getRoot(self):
		return self.root

	def insert(self,node,root):
		# Dynamically insert the node into the tree,
		# generating 'ghost' nodes along the way
		if root.id == node.id: #duplicate insertion
			if root.isGhost():
				root.ghost = False
			return self.updateWeights(root) #update weights and don't add node
		if root.id == "0":
			level = 0
		else:
			level = len(root.id)
		prefix = node.id[:level+1]
		for child in root.children:
			if child.id==prefix:
				return self.insert(node,child)
		#this subnode has no matching children,
		#do we need to create a ghost node or add this node?
		if prefix == node.id: #we're at the bottom of the tree, add the node
			root.children.append(node)
			node.parent = root
			self.nodes[node.id] = node
			self.edges[root.id+"-"+node.id] = 0
			return self.updateWeights(node)
		else: #we're not at the bottom, need a ghost node
			ghost = Node(prefix,root,True)
			root.children.append(ghost)
			self.nodes[ghost.id] = ghost
			self.edges[root.id+"-"+ghost.id] = 0
			return self.insert(node,ghost)

	def getNode(self, name):
		return self.nodes[name]

	def updateWeights(self,node):
		if node.parent == None:
			return
		self.edges[self.getEdge(node)] += 1
		return self.updateWeights(node.parent)

	def getEdge(self,node):
		return node.parent.id+"-"+node.id

	def export(self,treeFile):
		with open(treeFile, 'wb') as file:
			data = {"nodes":self.nodes,"edges":self.edges}
			pickle.dump(data, file)

	def loadFromFile(self,treeFile):
		with open(treeFile, 'rb') as file:
			data = pickle.load(file)
			self.nodes = data["nodes"]
			self.edges = data["edges"]