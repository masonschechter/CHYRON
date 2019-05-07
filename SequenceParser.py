import numpy as np
import openpyxl as xl
from pickle import load
from math import inf,log
class SequenceParser:

	def __init__(self, prefixMatch=0, magnitude = 100,decay_rate = 2):
		self.prefixMatch = prefixMatch
		self.magnitude = magnitude
		self.decay_rate = decay_rate
		# self.distFunc = self.EuclidianDist
		self.distFunc = self.linearDistance

	def string2path(self, seq):
		path = np.array([0.,0.,0.,0.]) # start at the origin
		m = self.magnitude
		for i, letter in enumerate(seq):
			# if not i % decay_change:
			# 	# change the decay rate
			# 	decay_rate 
			if letter == "A":
				path += np.array([1.,0.,0.,0.])*m
			elif letter == "C":
				path += np.array([0.,1.,0.,0.])*m
			elif letter == "G":
				path += np.array([0.,0.,1.,0.])*m
			elif letter == "T":
				path += np.array([0.,0.,0.,1.])*m
			# magnitude /= (1+ (1/(1+i)) )
			m /= self.decay_rate
			# if i % 3 == 0:
			# 	self.decay_rate += 1
		return path

	def stringDist(self,a,b):
		return self.distFunc(a,b)

	'''
		Euclidian Distance:
		Default function, takes the euclidian distance between two points
	'''
	def EuclidianDist(self,a,b):
		if a==b:
			return 0
		if not a[:self.prefixMatch] == b[:self.prefixMatch]:
			return -1
		p1 = self.string2path(a)
		p2 = self.string2path(b)
		return np.sqrt(np.sum((p1 - p2)**2))

	def linearDistance(self, a, b):
		if a==b:
			return 0 #Perfect score
		longerSeq, shorterSeq = (a, b) if len(a) >= len(b) else (b, a)
		highestScore = ( (1 - 2**(len(longerSeq)+1) ) / (-1) ) #Σ(4^n-1) over n=1 -> length of longer sequence
		prefix = 0
		for i in range(len(shorterSeq)):
			if not longerSeq[i] == shorterSeq[i]:
				break
			prefix += 1
		distance = ( 2 * prefix) / (len(longerSeq)+len(shorterSeq))
		#Okay now we have the prefix length.
		# if distance == 1:
		# 	return 0
		return 1-distance
		try:
			distance = 1 / ( log(1/(1-distance)) +0.01)
		except:
			print(f"distance: {distance} on seqs {a} and {b}")
			raise
		return distance

	def changeDistFunc(self, newFunc):
		self.distFunc = newFunc