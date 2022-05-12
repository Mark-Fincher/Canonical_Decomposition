"""
This class is intended to represent an undirected graph with edges labeled by positive integers.
We re-use the vertex and edge classes (vertex.py, edge.py) here. This might be a little confusing,
since these classes are made for quotient vertices and edges of simplicial orbifolds (or just 
triangulated manifolds). So of course many of their attributes and methods are meaningless 
when they're used here.

Here are the properties which matter.

example_vertex.Edges == List of edges having example_vertex as a vertex
example_edge.Vertices == List of the vertices of example_edge
example_edge.LocusOrder == The positive integer label for example_edge

The primary use of this data type is to represent the barycentric graph and singular locus
graph of a simplicial orbifold.
"""

class graph:
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