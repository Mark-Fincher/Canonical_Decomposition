"""
This file is used to create the barycentric graph and singular locus graph of an orbifold
triangulation. The graph data type is defined in graph.py.

The barycentric graph has one vertex for every tet, vertex class, edge class, and face class
of an orbifold triangulation "orb". You have an edge exactly when one of these is contained in the other.
There could be more than one edge between two vertices in the graph. The graph edges are given a
positive integer label corresponding to the order of the rotation group fixing it. The idea 
is that this graph is the 1-skeleton of the barycentric subdivision of the orbifold triangulation, 
quotiented out by the symmetries of the tets. 

For example, suppose 0-simplices V1 and V2 of tet0 belong to the same vertex class, V. If V1 
and V2 are not identified by a symmetry of tet0, then each contributes a distinct edge 
between the graph vertex of tet0 and the graph vertex of V in the barycentric graph. If they
are identified by a symmetry, then they only contribute one edge.

To say it more precisely, the orbit of V1 in tet0 under the action of the symmetry group of tet0
contributes exactly one edge from the graph vertex of tet0 to the graph vertex of V. You can make
the same kind of statement to determine graph edges between the other graph vertex types.

The singular locus graph is obtained from the barycentric graph by:
- Removing all edges labelled 1.
- Removing any isolated vertices which result from step 1.
- Removing any valence 2 vertex v, where the two connected edges are distint and having other vertex not equal to v.
  In this case, we can be sure that v is not an ideal point because of the classification of cusp types.

"""

from graph import*
from SimplicialOrbifold import*


# In addition to creating the barycentric graph, this function also sets CuspType for each
# vertex of orb. Must be one of the strings 'finite','torus','(2,2,2,2)','(2,3,6)','(2,4,4)',
# '(3,3,3)', or 'error'.
def barycentric_graph(orb):
	# Create the graph vertices and a dictionary associating them to
	# the tets, vertex classes, edge classes, and face classes.
	class_to_graph_vertex = dict()
	for tet in orb.Tetrahedra:
		class_to_graph_vertex[tet] = Vertex()
	for orb_vertex in orb.Vertices:
		class_to_graph_vertex[orb_vertex] = Vertex()
	for orb_edge in orb.Edges:
		class_to_graph_vertex[orb_edge] = Vertex()
	for face in orb.Faces:
		class_to_graph_vertex[face] = Vertex()
	graph_vertices = [class_to_graph_vertex[key] for key in class_to_graph_vertex.keys()]
	graph = Graph()
	graph.Vertices = graph_vertices
	# The rest of the function is about deciding what edges to make and what
	# their labels should be.
	# First we make the edges connected to the tet vertices.
	for tet in orb.Tetrahedra:
		vertex1 = class_to_graph_vertex[tet]
		skip = []
		for zero_subsimplex in ZeroSubsimplices:
			if zero_subsimplex not in skip:
				for sym in tet.Symmetries:
					skip.append(sym.image(zero_subsimplex))
				vertex2 = class_to_graph_vertex[tet.Class[zero_subsimplex]]
				if tet.non_trivial_sym_fixing(zero_subsimplex):
					graph.make_edge(vertex1,vertex2,3)
				else:
					graph.make_edge(vertex1,vertex2,1)
		skip = []
		for one_subsimplex in OneSubsimplices:
			if one_subsimplex not in skip:
				for sym in tet.Symmetries:
					skip.append(sym.image(one_subsimplex))
				vertex2 = class_to_graph_vertex[tet.Class[one_subsimplex]]
				if tet.non_trivial_sym_fixing(one_subsimplex):
					graph.make_edge(vertex1,vertex2,2)
				else:
					graph.make_edge(vertex1,vertex2,1)
		skip = []
		for two_subsimplex in TwoSubsimplices:
			if two_subsimplex not in skip:
				for sym in tet.Symmetries:
					skip.append(sym.image(two_subsimplex))
				vertex2 = class_to_graph_vertex[tet.Class[two_subsimplex]]
				if tet.non_trivial_sym_fixing(two_subsimplex):
					graph.make_edge(vertex1,vertex2,3)
				else:
					graph.make_edge(vertex1,vertex2,1)
	# Now we make edges connecting face vertices to edge and vertex vertices.
	for face in orb.Faces:
		vertex1 = class_to_graph_vertex[face]
		corner = face.get_glued_corner()
		tet = corner.Tetrahedron
		two_subsimplex = corner.Subsimplex
		# Record the face symmetries.
		face_syms = []
		for sym in tet.Symmetries:
			if sym.image(two_subsimplex) == two_subsimplex:
				face_syms.append(sym)
		if tet.face_glued_to_self(two_subsimplex):
			perm = tet.Gluing[two_subsimplex]
			extra_syms = [perm*sym for sym in face_syms]
			face_syms = face_syms + extra_syms
		skip = []
		for zero_subsimplex in ZeroSubsimplices:
			if is_subset(zero_subsimplex,two_subsimplex) and zero_subsimplex not in skip:
				for sym in face_syms:
					skip.append(sym.image(zero_subsimplex))
				vertex2 = class_to_graph_vertex[tet.Class[zero_subsimplex]]
				for sym in face_syms:
					if sym.image(zero_subsimplex) == zero_subsimplex and sym.tuple() != (0,1,2,3):
						graph.make_edge(vertex1,vertex2,2)
						break
				else:
					graph.make_edge(vertex1,vertex2,1)
		skip = []
		for one_subsimplex in OneSubsimplices:
			if is_subset(one_subsimplex,two_subsimplex) and one_subsimplex not in skip:
				for sym in face_syms:
					skip.append(sym.image(one_subsimplex))
				vertex2 = class_to_graph_vertex[tet.Class[one_subsimplex]]
				for sym in face_syms:
					if sym.image(one_subsimplex) == one_subsimplex and sym.tuple() != (0,1,2,3):
						graph.make_edge(vertex1,vertex2,2)
						break
				else:
					graph.make_edge(vertex1,vertex2,1)
	# Now connect edge vertices to vertex vertices.
	for edge in orb.Edges:
		vertex1 = class_to_graph_vertex[edge]
		# An edge either has a reflectional symmetry or no symmetry. It has the symmetry
		# if it's mapped nontrivially to itself by an order 2 sym of some tet or a face
		# being glued to itself. edge.has_symmetry() determines if this is the case.
		one_subsimplex = edge.Corners[0].Subsimplex
		tet = edge.Corners[0].Tetrahedron
		for zero_subsimplex in ZeroSubsimplices:
			if is_subset(zero_subsimplex,one_subsimplex):
				vertex2 = class_to_graph_vertex[tet.Class[zero_subsimplex]]
				graph.make_edge(vertex1,vertex2,edge.LocusOrder)
				if edge.has_symmetry():
					# In this case we want to stop the for loop, we only want to make
					# one new graph edge.
					break
	for i in range(len(graph.Vertices)):
		graph.Vertices[i].Index = i
	for i in range(len(graph.Edges)):
		graph.Edges[i].Index = i
	return graph

"""
Given an orbifold triang "orb", return its singular locus. First it gets the barycentric
graph, then it turns that into the singular locus.
"""
def singular_locus(orb):
	graph = barycentric_graph(orb)
	edges_list = [ edge for edge in graph.Edges ]
	for edge in edges_list:
		if edge.LocusOrder == 1:
			graph.delete_edge(edge)
	vertex_list = [ vertex for vertex in graph.Vertices]
	for vertex in vertex_list:
		if len(vertex.Edges) == 0:
			graph.delete_vertex(vertex)
		graph.attempt_remove_valence_2_vertex(vertex)
	return graph














				
