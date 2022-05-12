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

	def new_vertex(self):
		new = Vertex()
		new.Index = len(self.Vertices)
		self.Vertices.append(new)
		return new

	"""
	For vertices vertex1 and vertex2 already in the graph, create a new edge connecting 
	vertex1 to vertex 2 with integer label "label".
	"""
	def make_edge(self,vertex1,vertex2,label):
		newEdge = Edge()
		newEdge.Index = len(self.Edges)
		self.Edges.append(newEdge)
		newEdge.Vertices = [vertex1,vertex2]
		newEdge.LocusOrder = label
		vertex1.Edges.append(newEdge)
		if vertex1 != vertex2:
			vertex2.Edges.append(newEdge)

	"""
	Delete an edge from the graph.  
	"""
	def delete_edge(self,edge):
		vertex0 = edge.Vertices[0]
		vertex1 = edge.Vertices[1]
		vertex0.Edges.remove(edge)
		if edge in vertex1.Edges:
			vertex1.Edges.remove(edge)
		self.Edges.remove(edge)
		edge.Vertices = None
		for i in range(len(self.Edges)):
			self.Edges[i].Index = i

	"""
	Delete a vertex from the graph. This requires that we delete all edges containing
	that vertex.
	"""
	def delete_vertex(self,vertex):
		edges_to_delete = [ edge for edge in vertex.Edges]
		for edge in edges_to_delete:
			self.delete_edge(edge)
		self.Vertices.remove(vertex)
		for i in range(len(self.Vertices)):
			self.Vertices[i].Index = i

	"""
	IF a vertex has valence 2, by which I mean it has exactly two edges connected to it,
	both with other endpoints connected to other vertices than this one. AND IF both those edges
	have the same label. THEN we make the two edges into one edge, removing the vertex.

	This function sees if a given vertex satisfies the above criteria, and if so then it 
	does the modification. In that case it returns True, otherwise False.
	"""
	def attempt_remove_valence_2_vertex(self,vertex):
		if len(vertex.Edges) != 2:
			return False
		other_vertex = []
		for edge in vertex.Edges:
			vertex0 = edge.Vertices[0]
			vertex1 = edge.Vertices[1]
			if vertex1 != vertex:
				other_vertex.append(vertex1)
			elif vertex0 != vertex:
				other_vertex.append(vertex0)
			else:
				return False
		edge0 = vertex.Edges[0]
		edge1 = vertex.Edges[1]
		if edge0.LocusOrder != edge1.LocusOrder:
			return False
		label = edge0.LocusOrder
		self.delete_vertex(vertex)
		self.make_edge(other_vertex[0],other_vertex[1],label)
		return True


