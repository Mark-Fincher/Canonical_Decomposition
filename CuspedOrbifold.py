"""
Cusped orbifold class. Like the t3m triangulation class, except we only include stuff we need,
and our triangulations are of orbifolds. Much of this is copied from Goerner et al, and that is
built on snappy.
"""

from HoroTriangle import*
from simplex import *
from tetrahedron import Tetrahedron
from corner import Corner
from arrow import Arrow
#from .face import Face
from edge import Edge
from vertex import Vertex
from Exact_Arithmetic import*

class CuspedOrbifold:
	def __init__(self,tets_list):
		self.Tetrahedra = tets_list
		self.Edges = []
		self.Faces = []
		self.Vertices = []
		self.is_canonical = None
		self.DestSeq = None
		self.PachnerPath = []
		for T in self.Tetrahedra:
			T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
		self.build_vertex_classes()
		self.build_edge_classes()
		self.add_cusp_cross_sections()
		for cusp in self.Vertices:
			self.normalize_cusp(cusp)
		self.see_if_canonical()

	def add_cusp_cross_sections(self):
		for cusp in self.Vertices:
			self.add_one_cusp_cross_section(cusp)

	def add_one_cusp_cross_section(self, cusp):
		"Build a cusp cross section as described in Section 3.6 of the paper"
		tet0, vert0 = tets_and_vertices_of_cusp(cusp)[0]
		face0 = FacesAnticlockwiseAroundVertices[vert0][0]
		tet0.horotriangles[vert0] = HoroTriangle(
			tet0, vert0, face0,
			SquareRootCombination.One())
		active = [(tet0, vert0)]
		while active:
			tet0, vert0 = active.pop()
			for face0 in FacesAnticlockwiseAroundVertices[vert0]:
				# Mark added the following line (and only that line) 8/2/2021. Because some faces might be unglued for us.
				if tet0.Neighbor[face0] != None:
					tet1, face1 = glued_to(tet0, face0)
					vert1 = tet0.Gluing[face0].image(vert0)
					if tet1.horotriangles[vert1] is None:
						tet1.horotriangles[vert1] = HoroTriangle(tet1, vert1, face1,
								tet0.horotriangles[vert0].lengths[face0])
						active.append( (tet1, vert1) )
			# Following added by Mark 8/9/2021 to acommodate symmetries in building horotriangles
			for perm in tet0.Symmetries:
				if perm.image(vert0) != vert0:
					vert1 = perm.image(vert0)
					if tet0.horotriangles[vert1] is None:
						face_0 = FacesAnticlockwiseAroundVertices[vert0][0]
						face1 = perm.image(face_0)
						tet0.horotriangles[vert1] = HoroTriangle(tet0,vert1,face1,tet0.horotriangles[vert0].lengths[face_0])
						active.append((tet0,vert1))



	def _get_cusp(self, cusp):
		"""
		Helper method so the user can specify a cusp by its index as well
		the actual t3m.Vertex.
		"""
		if not isinstance(cusp, Vertex):
			cusp = self.Vertices[cusp]
		return cusp


	"""
	Following changed by Mark, 7/12/2021. In the orbifold case, the cusp area computation is a little different.
	If a horotriangle is at a vertex which is fixed by a nontrivial symmetry, then the area of the triangle is
	divided by 3. And if a vertex is mapped to another vertex of the same tetrahedron by a symmetry, then only
	one of the corresponding horotriangles contributes area to the cusp.
	"""
	def cusp_area(self, cusp):
		cusp = self._get_cusp(cusp)
		area = SquareRootCombination.Zero()
		already_counted = []
		for T, V in tets_and_vertices_of_cusp(cusp):
			if (T,V) not in already_counted:
				already_counted.append((T,V))
				V_index = ZeroSubsimplices.index(V)
				fix_V = []
				for perm in T.Symmetries:
					if perm[V_index] == V_index:
						fix_V.append(perm)
					if (T,ZeroSubsimplices[perm[V_index]]) not in already_counted:
						already_counted.append((T,ZeroSubsimplices[perm[V_index]]))
				if len(fix_V) > 1:
					area += T.horotriangles[V].area/SquareRootCombination([(1,3)])
				else:
					area += T.horotriangles[V].area
		return area

	def rescale_cusp(self, cusp, scale):
		cusp = self._get_cusp(cusp)
		for T, V in tets_and_vertices_of_cusp(cusp):
			T.horotriangles[V].rescale(scale)


	def normalize_cusp(self, cusp):
		"""
		Rescale cusp to have area sqrt(3). This choice ensures that
		all tilts are again Q-linear combinations of square roots
		of integers.
		"""
		cusp = self._get_cusp(cusp)

		target_area = SquareRootCombination.SqrtThree()

		area = self.cusp_area(cusp)
		ratio = (target_area/area).sqrt()
		self.rescale_cusp(cusp, ratio)

	def LHS_of_convexity_equations(self):
		"""
		For each face in the triangulation, return a quantity which is < 0
		if and only if the corresponding pair of tetrahedra are strictly
		convex.
		"""
		ans = []
		for tet0 in self.Tetrahedra:
			for vert0 in ZeroSubsimplices:
				tet1, face1 = glued_to(tet0, comp(vert0))
				ans.append(tet0.tilt(vert0) + tet1.tilt(comp(face1)))
		return ans

	def find_opaque_faces(self):
		"""
		Returns a list of bools indicating whether a face of a tetrahedron
		of the given proto-canonical triangulation is opaque.
		The list is of the form
		[ face0_tet0, face1_tet0, face2_tet0, face3_tet0, face0_tet1, ...]
		"""
		num_cusps = len(self.Vertices)
		for i in range(num_cusps):
			self.normalize_cusp(i)
        
		tilts = self.LHS_of_convexity_equations()
		result = []
		for tilt in tilts:
			# Face is transparent when tilt is exactly 0
			if tilt == SquareRootCombination.Zero():
				result.append(False)
			# Use interval aritmetic to certify tilt is negative
			elif tilt.evaluate() < 0:
				result.append(True)
			else:
				# We failed
				raise Exception(
					"Could not certify proto-canonical triangulation")
                
		return result

    

    # Construct the vertices.
	#
	def build_vertex_classes(self):
		for tet in self.Tetrahedra:
			for zero_subsimplex in ZeroSubsimplices:
				if ( tet.Class[zero_subsimplex] == None ):
					newVertex = Vertex()
					self.Vertices.append(newVertex)
					self.walk_vertex(newVertex,zero_subsimplex,tet)
		for i in range(len(self.Vertices)):
			self.Vertices[i].Index = i

	def walk_vertex(self,vertex,zero_subsimplex,tet):
		if (tet.Class[zero_subsimplex] != None ):
			return
		else:
			tet.Class[zero_subsimplex] = vertex
			vertex.Corners.append(Corner(tet,zero_subsimplex))
			for two_subsimplex in TwoSubsimplices:
				if ( is_subset(zero_subsimplex,two_subsimplex)
                     and
                     tet.Gluing[two_subsimplex] != None):
						self.walk_vertex(vertex,
                                     tet.Gluing[two_subsimplex].image(zero_subsimplex),
                                     tet.Neighbor[two_subsimplex])
			# Following added by Mark 8/9/2021, if two vertices in a tet are identified by a symmetry of
			# the tet, they should be in the same vertex class.
			for perm in tet.Symmetries:
				if perm.image(zero_subsimplex) != zero_subsimplex:
					self.walk_vertex(vertex,perm.image(zero_subsimplex),tet)

	
	"""
	Two edges are in the same edge class if they're identified by a face pairing or symmetry.
	We want to determine all the edge classes. We also want to determine if an edge is part of the
	singular locus or not. And if it is, what is the order of that local group? This function figures
	out all of these things.

	Getting the local group order is a little tricky. You can move around a quotient edge in multiple
	ways, because of the symmetries of some tets. If one path takes you fully around a quotient edge
	in pi/3 radians, and another path takes you around in pi/6 radians, then the group element of the
	first is the square of the group element of the second. Need to find the path with smallest total
	angle; that will correspond to the generator of the local group.

	Each quotient edge corresponds to a distinct Edge object. The corners of the Edge object
	correspond to a particular shortest walk around the quotient edge. So not every one-simplex
	which quotients to that edge will be in the corners list. You can see which Edge object
	a one-simplex is assigned as Tet.Class[one_simplex].
	"""
	def build_edge_classes(self):
		for tet in self.Tetrahedra:
			for one_subsimplex in OneSubsimplices:
				if tet.Class[one_subsimplex] is None:
					a = Arrow(one_subsimplex,RightFace[one_subsimplex],tet)
					first_arrow = a.copy()
					completed_walks = []
					active = [[a]]
					while active:
						walk = active.pop()
						a = walk[-1]
						for sym in a.Tetrahedron.Symmetries:
							new_a = Arrow(sym.image(a.Edge),sym.image(a.Face),a.Tetrahedron)
							new_a.next()
							# if new_a is None, we don't do anything.
							if new_a != None:
								# now we want to check if new_a, or the image under some symmetry
								# of new_a, is first_arrow. If so, then walk is a complete walk.
								complete = False
								if new_a.Tetrahedron == first_arrow.Tetrahedron:
									for other_sym in new_a.Tetrahedron.Symmetries:
										if other_sym.image(new_a.Edge) == first_arrow.Edge and other_sym.image(new_a.Face) == first_arrow.Face:
											complete = True
								if complete:
									completed_walks.append(walk)
								else:
									walk.append(new_a)
									active.append(walk)
					shortest_walk = completed_walks[0]
					for walk in completed_walks:
						if len(walk) < len(shortest_walk):
							shortest_walk = walk
					new_edge = Edge()
					z = ComplexSquareRootCombination.One()
					for a in shortest_walk:
						z = z*a.Tetrahedron.edge_params[a.Edge]
						new_edge._add_corner(a)
						for sym in a.Tetrahedron.Symmetries:
							a.Tetrahedron.Class[sym.image(a.Edge)] = new_edge
					# Now find n s.t. z^n = 1
					n = 0
					w = ComplexSquareRootCombination.One()
					while n < 1000:
						w = w*z
						n = n + 1
						# So w = z^n
						if w == ComplexSquareRootCombination.One():
							break
						if n == 999:
							print("error in LocusOrder")
					new_edge.LocusOrder = n
					self.Edges.append(new_edge)
		for i in range(len(self.Edges)):
			self.Edges[i].Index = i
					








	"""
	def build_edge_classes(self):
		for tet in self.Tetrahedra:
			for one_subsimplex in OneSubsimplices:
				if tet.Class[one_subsimplex] == None:
					newEdge = Edge()
					self.Edges.append(newEdge)
					first_arrow = Arrow(one_subsimplex,RightFace[one_subsimplex],tet)
					a = first_arrow.copy()
					sanity_check = 0
					while 1:
						print(a)
						if sanity_check > 6*len(self.Tetrahedra):
							print("error building edge classes")
							break
						newEdge._add_corner(a)
						a.Tetrahedron.Class[a.Edge] = newEdge
						if a.next() == None:
							print("starting from edge",SubsimplexName[one_subsimplex])
							print("the following arrow has no next",a)
							# For us, that doesn't mean we hit the boundary, because there is
							# no boundary. It just means we need to apply a symmetry to get the
							# face gluing data. In this case, the arrow a was not changed, though
							# typically a.next() does change a.
							for sym in a.Tetrahedron.Symmetries:
								if a.Tetrahedron.Neighbor[sym.image(a.Face)] != None:
									print(sym)
									a = Arrow(sym.image(a.Edge),sym.image(a.Face),a.Tetrahedron)
									print(a)
									a.next()
									print(a)
									break
						if a == first_arrow:
							break
						sanity_check = sanity_check + 1
		for i in range(len(self.Edges)):
			self.Edges[i].Index = i
		for edge in self.Edges:
			z = ComplexSquareRootCombination.One()
			for corner in edge.Corners:
				z = z*corner.Tetrahedron.edge_params[corner.Subsimplex]
			# Now find n s.t. z^n = 1
			n = 0
			w = ComplexSquareRootCombination.One()
			while n < 1000:
				w = w*z
				n = n + 1
				# So w = z^n
				if w == ComplexSquareRootCombination.One():
					break
				if n == 999:
					print("error in LocusOrder")
			edge.LocusOrder = n
	"""







	# set self.is_canonical to True or False
	def see_if_canonical(self):
		for tet1 in self.Tetrahedra:
			for face1 in TwoSubsimplices:
				if tet1.Neighbor[face1] != None:
					tet2 = tet1.Neighbor[face1]
					face2 = tet1.Gluing[face1].image(face1)
					if (tet1.tilt(comp(face2)) + tet2.tilt(comp(face2))).evaluate() > 0:
						self.is_canonical = False
						return
		# if it hasn't already returned, set it to True.
		self.is_canonical = True



    