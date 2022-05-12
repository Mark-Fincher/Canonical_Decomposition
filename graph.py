"""
This class is intended to represent an undirected graph with edges labeled by positive integers.
We re-use the vertex and edge classes (vertex.py, edge.py) here. This might be a little confusing,
since these classes are made for quotient vertices and edges of simplicial orbifolds (or just 
triangulated manifolds). So of course many of their attributes and methods are meaningless 
when they're used here.

Here are the properties which matter.

example_vertex.Edges == List of the distinct edges having example_vertex as a vertex
example_edge.Vertices == List of the (not necessarily distinct) vertices of example_edge
example_edge.LocusOrder == The positive integer label for example_edge

The primary use of this data type is to represent the barycentric graph and singular locus
graph of a simplicial orbifold.
"""

from vertex import*
from edge import*

class Graph:
	def __init__(self):
		self.Vertices = []
		self.Edges = []

	def info(self):
		print('Vertices = ',self.Vertices)
		print(' ')
		print('Edges = ',self.Edges)
		print(' ')
		for edge in self.Edges:
			print('Label of edge',edge,'is',edge.LocusOrder)
			print('Vertices of edge',edge,'are',edge.Vertices)
			print(' ')

	# For vertices vertex1 and vertex2 already in the graph, create a new edge connecting 
	# vertex1 to vertex 2 with integer label "label".
	def make_edge(self,vertex1,vertex2,label):
		newEdge = Edge()
		self.Edges.append(newEdge)
		newEdge.Vertices = [vertex1,vertex2]
		newEdge.LocusOrder = label
		vertex1.Edges.append(newEdge)
		if vertex1 != vertex2:
			vertex2.Edges.append(newEdge)
