"""
Cusped orbifold class. Similar to the SnapPy t3m mcomplex class. Some things are taken, perhaps
with modification, from ancillary files of the paper "A census of tetrahedral hyperbolic 
manifolds" by Evgeny Fominykh, Stavros Garoufalidis, Matthias Goerner, Vladimir Tarkaev, 
and Andrei Vesnin. I comment FGGTV before those things.
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
		# Make the indices of the tetrahedra match their indices in the list.
		for i in range(len(self.Tetrahedra)):
			self.Tetrahedra[i].Index = i
		# Make sure all tets have at least the identity map in their Symmetries list, required by some functions.
		for tet in self.Tetrahedra:
			if len(tet.Symmetries) == 0:
				tet.Symmetries.append(Perm4((0,1,2,3)))
		# Self is geometric if it has values assigned to its edge parameters. 
		if self.Tetrahedra[0].edge_params[E01] is None:
			self.complex_arithmetic_type = None
			self.real_arithmetic_type = None
			self.is_geometric = False
		else:
			self.complex_arithmetic_type = type(self.Tetrahedra[0].edge_params[E01])
			self.real_arithmetic_type = type(self.Tetrahedra[0].edge_params[E01].real)
			self.is_geometric = True
		self.Edges = []
		self.Faces = []
		self.Vertices = []
		self.DestSeq = None
		self.PachnerPath = []
		for T in self.Tetrahedra:
			T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
		self.build_vertex_classes()
		self.build_edge_classes()
		if self.is_geometric:
			self.add_cusp_cross_sections()
			for cusp in self.Vertices:
				self.normalize_cusp(cusp)

	def info(self):
		tets = self.Tetrahedra
		if self.is_geometric:
			print('this triangulation is geometric')
		else:
			print('this triangulation is not geometric')
		print('number of tetrahedra is',len(tets))
		for i in range(len(tets)):
			print('gluing data for', tets[i], 'is')
			for j in range(4):
				print('face', j, 'of', tets[i], 'glued to', tets[i].Neighbor[TwoSubsimplices[j]])
				print('with gluing map',tets[i].Gluing[TwoSubsimplices[j]])
			print('symmetries of',tets[i],'are')
			for sym in tets[i].Symmetries:
				print(sym)
			if self.is_geometric:
				print('edge_params of',tets[i],'are')
				print('edge 01',tets[i].edge_params[E01].real,'+',tets[i].edge_params[E01].imag,'* i')
				print('edge 02',tets[i].edge_params[E02].real,'+',tets[i].edge_params[E02].imag,'* i')
				print('edge 03',tets[i].edge_params[E03].real,'+',tets[i].edge_params[E03].imag,'* i')
				print('edge 12',tets[i].edge_params[E12].real,'+',tets[i].edge_params[E12].imag,'* i')
				print('edge 13',tets[i].edge_params[E13].real,'+',tets[i].edge_params[E13].imag,'* i')
				print('edge 23',tets[i].edge_params[E23].real,'+',tets[i].edge_params[E23].imag,'* i')
			print('edge labels of',tets[i],'are')
			print('edge 01: ',tets[i].edge_labels[E01])
			print('edge 02: ',tets[i].edge_labels[E02])
			print('edge 03: ',tets[i].edge_labels[E03])
			print('edge 12: ',tets[i].edge_labels[E12])
			print('edge 13: ',tets[i].edge_labels[E13])
			print('edge 23: ',tets[i].edge_labels[E23])

	def copy(self):
		new_tets = []
		new_to_old = {}
		old_to_new = {}
		for tet in self.Tetrahedra:
			new_tet = Tetrahedron()
			old_to_new[tet] = new_tet
			new_to_old[new_tet] = tet
			new_tets.append(new_tet)
		for new_tet in new_tets:
			new_tet.Index = new_to_old[new_tet].Index
			for one_subsimplex in OneSubsimplices:
				if self.is_geometric:
					new_tet.edge_params[one_subsimplex] = new_to_old[new_tet].edge_params[one_subsimplex]
				new_tet.edge_labels[one_subsimplex] = new_to_old[new_tet].edge_labels[one_subsimplex]
			for sym in new_to_old[new_tet].Symmetries:
				new_tet.Symmetries.append(sym)
			for face in TwoSubsimplices:
				if new_to_old[new_tet].Neighbor[face] != None:
					new_tet.attach(face,old_to_new[new_to_old[new_tet].Neighbor[face]],
						new_to_old[new_tet].Gluing[face].tuple())
		return CuspedOrbifold(new_tets)

	def remove_geometry(self):
		"""
		If self is geometric, remove its geometry by setting edge_params and horotriangles to None and is_geometric
		to False. If self is already not geometric, do nothing.
		"""
		if self.is_geometric is False:
			return
		for tet in self.Tetrahedra:
			tet.edge_params = {E01:None,E23:None,E02:None,E13:None,E03:None,E12:None}
        	tet.horotriangles = {V0:None, V1:None, V2:None, V3:None}
        self.is_geometric = False

	"""
	Clear all classes and horotriangles, then rebuild this data. Want to do this after changing
	the triangulation.
	"""
	def clear_and_rebuild(self):
		self.Edges = []
		self.Vertices = []
		for tet in self.Tetrahedra:
			tet.clear_Class()
			tet.horotriangles = {V0:None, V1:None, V2:None, V3:None}
		for i in range(len(self.Tetrahedra)):
			self.Tetrahedra[i].Index = i
		self.build_vertex_classes()
		self.build_edge_classes()
		self.add_cusp_cross_sections()
		for cusp in self.Vertices:
			self.normalize_cusp(cusp)

	"""
	The following is a very important method which is used at the end of every Pachner move. It updates
	the edge and vertex classes, horotriangles, and tet indices.

	We do this in a very simple way for the indices. Just label the tets according to their indices in self.Tetrahedra.

	For edge classes and horotriangles, it is also simple. We just throw out the previous data
	we had for them and re-build them from scratch (assuming, for horotriangles, we already have updated vertex class data).

	For vertex classes, we choose to make things more complicated. We could just throw them out and rebuild with
	self.build_vertex_classes(), but there are reasons to try to instead update existing vertex classes. I.e. you might
	want to refer to some particular vertex class (cusp) before and after a series of Pachner moves. Therefore, we keep
	the same vertex classes and just update their corners. Meaning if a tet was removed, its 0-simplices should be removed
	the corners lists of whatever vertex classes they belonged to. And if a new tet was added in, its 0-simplices should
	be added to the correct vertex corners lists and its Class dictionary should be accordingly updated.

	An important point: in the case of a CuspedOrbifold Pachner move, vertex classes cannot be created nor destroyed (because
	vertex classes correspond exactly to cusps). But we DO have to worry about that for a non-geometric Pachner move (i.e. for
	SimplicialOrbifold objects).
	"""
	def clear_and_rebuild_new(self):
		# Reset the tet indices.
		for i in range(len(self.Tetrahedra)):
			self.Tetrahedra[i].Index = i
		# Clear edge classes.
		self.Edges = []
		for tet in self.Tetrahedra:
			tet.horotriangles = {V0:None, V1:None, V2:None, V3:None}
			for one_subsimplex in OneSubsimplices:
				tet.Class[one_subsimplex] = None
		# Rebuild edge classes.
		self.build_edge_classes()
		# Now we try to make the vertex classes correct. Because of what we did in the Pachner move, it should currently be the case that
		# for every tet and every zero_subsimplex, tet.Class[zero_subsimplex] is the correct vertex class. The only potentially incorrect thing
		# at this moment is that each vertex class could have the wrong set of corners. So we need to fix that.
		Tets_ZeroSubsimplices = []
		for tet in self.Tetrahedra:
			for zero_subsimplex in ZeroSubsimplices:
				Tets_ZeroSubsimplices.append((tet,zero_subsimplex))
		for vertex in self.Vertices:
			vertex.Corners = [Corner(pair[0],pair[1]) for pair in Tets_ZeroSubsimplices if pair[0].Class[pair[1]] is vertex]
		# Let's do a check that the vertex classes are correct.
		self.check_vertex_classes()
		# Now re-build horotriangles/cusp geometry. Note that we couldn't do this until we had correct vertex classes.
		if self.is_geometric:
			self.add_cusp_cross_sections()
			for cusp in self.Vertices:
				self.normalize_cusp(cusp)

	# Check the vertex classes are correctly assigned.
	def check_vertex_classes(self):
		seen_vertex_classes = []
		for tet in self.Tetrahedra:
			for zero_subsimplex in ZeroSubsimplices:
				seen_vertex_class.append(tet.Class[zero_subsimplex])
				for sym in tet.Symmetries:
					assert tet.Class[zero_subsimplex] == tet.Class[sym.image(zero_subsimplex)]
				for two_subsimplex in TwoSubsimplices:
					if is_subset(zero_subsimplex,two_subsimplex) and tet.Neighbor[two_subsimplex] is not None:
						nbr = tet.Neighbor[two_subsimplex]
						gluing = tet.Gluing[two_subsimplex]
						assert tet.Class[zero_subsimplex] == nbr.Class[gluing.image(zero_subsimplex)]
		assert set(seen_vertex_class) == set(self.Vertices)

	def add_tet(self, tet):
		self.Tetrahedra.append(tet)

	# Add one new tetrahedron and return one of its arrows.
	def new_arrow(self):
		tet = Tetrahedron()
		self.add_tet(tet)
		return Arrow(E01,F3,tet)

	# Or, add a whole bunch of them.
	#
	def new_arrows(self,n):
		return [self.new_arrow() for i in range(n)]

	# FGGTV
	def add_cusp_cross_sections(self):
		for cusp in self.Vertices:
			self.add_one_cusp_cross_section(cusp)

	# FGGTV
	def add_one_cusp_cross_section(self, cusp):
		"Build a cusp cross section as described in Section 3.6 of the paper"
		tet0, vert0 = tets_and_vertices_of_cusp(cusp)[0]
		face0 = FacesAnticlockwiseAroundVertices[vert0][0]
		tet0.horotriangles[vert0] = HoroTriangle(
			tet0, vert0, face0,
			self.real_arithmetic_type.One())
		active = [(tet0, vert0)]
		while active:
			tet0, vert0 = active.pop()
			for face0 in FacesAnticlockwiseAroundVertices[vert0]:
				# We need the following for orbifolds because some faces might be unglued.
				if tet0.Neighbor[face0] != None:
					tet1, face1 = glued_to(tet0, face0)
					vert1 = tet0.Gluing[face0].image(vert0)
					if tet1.horotriangles[vert1] is None:
						tet1.horotriangles[vert1] = HoroTriangle(tet1, vert1, face1,
								tet0.horotriangles[vert0].lengths[face0])
						active.append( (tet1, vert1) )
			# With orbifolds, need to acommodate symmetries in building horotriangles.
			for perm in tet0.Symmetries:
				if perm.image(vert0) != vert0:
					vert1 = perm.image(vert0)
					if tet0.horotriangles[vert1] is None:
						face_0 = FacesAnticlockwiseAroundVertices[vert0][0]
						face1 = perm.image(face_0)
						tet0.horotriangles[vert1] = HoroTriangle(tet0,vert1,face1,tet0.horotriangles[vert0].lengths[face_0])
						active.append((tet0,vert1))



	# FGGTV
	def _get_cusp(self, cusp):
		"""
		Helper method so the user can specify a cusp by its index as well
		the actual t3m.Vertex.
		"""
		if not isinstance(cusp, Vertex):
			cusp = self.Vertices[cusp]
		return cusp


	"""
	FGGTV. In the orbifold case, the cusp area computation is a little different.
	If a horotriangle is at a vertex which is fixed by a nontrivial symmetry, then 
	the area of the triangle is divided by 3. And if a vertex is mapped to another 
	vertex of the same tetrahedron by a symmetry, then only one of the corresponding 
	horotriangles contributes area to the cusp.
	"""
	def cusp_area(self, cusp):
		cusp = self._get_cusp(cusp)
		area = self.real_arithmetic_type.Zero()
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

	# FGGTV
	def rescale_cusp(self, cusp, scale):
		cusp = self._get_cusp(cusp)
		for T, V in tets_and_vertices_of_cusp(cusp):
			T.horotriangles[V].rescale(scale)


	# FGGTV
	def normalize_cusp(self, cusp):
		"""
		Rescale cusp to have area sqrt(3). This choice ensures that
		all tilts are again Q-linear combinations of square roots
		of integers.
		"""
		cusp = self._get_cusp(cusp)

		One = self.real_arithmetic_type.One()
		Three = One + One + One
		# This target area is chosen because the arithmetic works nicely for orbifolds commensurable
		# to the fig eight knot complement. See FGGTV. It works fine for other orbifolds too if
		# we are just using sage interval arithmetic.
		target_area = Three.sqrt()

		area = self.cusp_area(cusp)
		ratio = (target_area/area).sqrt()
		self.rescale_cusp(cusp, ratio)

	# FGGTV
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

	# FGGTV
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
			# ORBIFOLDS. If two vertices in a tet are identified by a symmetry of
			# the tet, they should be in the same vertex class.
			for perm in tet.Symmetries:
				if perm.image(zero_subsimplex) != zero_subsimplex:
					self.walk_vertex(vertex,perm.image(zero_subsimplex),tet)

	"""
	Two edges are in the same class if they're identified by a sequence of face gluings and
	symmetries. In that case, the two edges are assigned the same edge object (defined in edge.py).
	For example, if edge E01 of tet0 has the edge object EDGE, then tet0.Class[E01] == EDGE.
	One attribute of the edge class is "Corners". A Corner is an object with attributes Tetrahedron
	and Subsimplex. It might make sense for EDGE.Corners to be the list of all Corners such that
	Corner.Tetrahedron.Class[Corner.Subsimplex] == EDGE. But in fact we do things slightly
	differently. We can walk around EDGE using arrows, and only the (one)Subsimplex-Tetrahedron pairs we
	encounter in a complete walk get recorded into EDGE.Corners. Because of symmetries, this need not be all of the
	1-simplices belonging to this edge class. Additionally, there could be repeats in this list, if an 
	arrow and its reverse both apppear in the walk. 

	We do it this way because, in the geometric case, we need to know the total dihedral angle traversed in a
	walk, theta. Then the integer n satisfying n*theta = 2*pi is the order of the local group fixing the
	edge, recorded as EDGE.LocusOrder == n.

	If self is not geometric, then of course we don't worry about this. We just set EDGE.LocusOrder to be
	the edge_label of any of the 1-simplices belonging to EDGE.

	From this perpective, it would really make more sense for EDGE.Corners to be a list of arrows in
	a complete walk around EDGE. We choose not to do this, but that's the way to think of it. 
	"""
	def build_edge_classes(self):
		for tet in self.Tetrahedra:
			for one_subsimplex in OneSubsimplices:
				if tet.Class[one_subsimplex] is None:
					newEdge = Edge()
					self.Edges.append(newEdge)
					a = Arrow(one_subsimplex,RightFace[one_subsimplex],tet)
					first_arrow = a.copy()
					sanity_check = 0
					unfinished_walk = True
					while unfinished_walk:
						if sanity_check > 12*len(self.Tetrahedra):
							raise Exception('Bad gluing data: could not construct edge link.')
						newEdge.Corners.append(Corner(a.Tetrahedron, a.Edge))
						if a.true_next() is None:
							raise Exception('Hit boundary. Did not construct edge link.')
						else:
							#Check if a is now first_arrow, or the image of a under 
							#a symmetry is first_arrow.
							if a.Tetrahedron is first_arrow.Tetrahedron:
								for sym in a.Tetrahedron.Symmetries:
									if ( sym.image(a.Edge) == first_arrow.Edge 
										and 
										sym.image(a.Face) == first_arrow.Face):
										#Then we've finished walking around the edge,
										#and we want to stop the while loop.
										unfinished_walk = False
										break
							sanity_check = sanity_check + 1
					#We've finished walking around the edge, and now we want to assign newEdge to all
					#edges in this class. This might not just be the edges in newEdge.Corners; have
					#to apply symmetries to get other edges in the class.
					for corner in newEdge.Corners:
						for sym in corner.Tetrahedron.Symmetries:
							corner.Tetrahedron.Class[sym.image(corner.Subsimplex)] = newEdge
					#Now we figure out the order of the local group fixing newEdge.
					if newEdge.Corners[0].Tetrahedron.edge_labels[Corners.Subsimplex] is not None:
						# Then just use this edge label to set the LocusOrder.
						newEdge.LocusOrder = newEdge.Corners[0].Tetrahedron.edge_labels[Corners.Subsimplex]
					elif self.is_geometric:
						z = self.complex_arithmetic_type.One()
						for corner in newEdge.Corners:
							z = z*corner.Tetrahedron.edge_params[corner.Subsimplex]
						# LocusOrder is the smallest positive integer n such that z^n == 1.
						n = 0
						w = self.complex_arithmetic_type.One()
						while n < 1000:
							w = w*z
							n = n + 1
							# So w = z^n
							if w == self.complex_arithmetic_type.One():
								newEdge.LocusOrder = n
								break
						else:
							raise Exception('error getting edge locus order')
					else:
						raise Exception('cannot determine LocusOrder of edge because the orb is not geometric and has no edge_labels')
		# Make sure all edge labels are assigned. If we initiated a CuspedOrbifold object with edge parameters but not edge labels, then
		# edge labels would not be assigned at this point.
		for tet in self.Tetrahedra:
			for one_subsimplex in OneSubsimplices:
				if tet.edge_labels[one_subsimplex] is None:
					tet.edge_labels[one_subsimplex] = tet.Class[one_subsimplex].LocusOrder
		# Set the indices of the edges.				
		for i in range(len(self.Edges)):
			self.Edges[i].Index = i

	# Checks if the triangulation is proto_canonical or not.
	def is_proto_canonical(self):
		# Can't talk about being proto-canonical if it's not geometric.
		assert self.is_geometric
		for tet1 in self.Tetrahedra:
			if tet1.is_flat():
				#don't try to compute tilts of a flat tetrahedron.
				continue
			for face1 in TwoSubsimplices:
				if tet1.Neighbor[face1] is not None:
					tet2 = tet1.Neighbor[face1]
					if tet2.is_flat():
						#again, don't try to compute tilts of a flat tetrahedron.
						continue
					face2 = tet1.Gluing[face1].image(face1)
					if (tet1.tilt(comp(face1)) + tet2.tilt(comp(face2))).evaluate() > 0:
						return False
		# Now make sure that if there are any flat tetrahedra, they are admissible.
		for tet in self.Tetrahedra:
			if tet.is_flat() and self.check_admissible(tet) == False:
				return False
		return True

	"""
	A flat tetrahedron is "admissible" if it represents a square face being glued to itself,
	and if doing a special 4-4 move there results in a concave edge. See the write-up for a
	more detailed explanation. The point is that our orbifold proto-canonical triangulations
	very well could have flat tets in them (which I think is not the case for manifolds).
	But you can't compute the tilt of a flat tet, so how are you really sure it's proto-canonical?
	This is my solution. If a flat tet is "admissible", then it belongs. Otherwise the triangulation
	is not proto-canonical.

	Should update this. Just have it compute tilts instead of doing special_four_to_four.
	"""
	def check_admissible(self,tet):
		if tet.is_flat() is False:
			raise Exception("Tried to to check if a non-flat tet is admissible.")
		tet_index = tet.Index
		copy = self.copy()
		tet = copy.Tetrahedra[tet_index]
		for one_subsimplex in OneSubsimplices:
			if copy.special_four_to_four(tet.Class[one_subsimplex]):
				# Because of how special_four_to_four is written, the three new tets which
				# make up the octahedron are in positions -1, -2, and -3 of copy.Tetrahedra.
				# We know the interior edge of the newly created octahedron will be E01 of any
				# of those tetrahedra. This tet is admissible iff this edge class is concave.
				# To check it it's concave, we check the concavity of the faces adjacent to it.
				# They are F2 and F3 of copy.Tetrahedra[-2].
				new_0 = copy.Tetrahedra[-3]
				new_1 = copy.Tetrahedra[-2]
				new_2 = copy.Tetrahedra[-1]
				if new_0.is_flat() and new_2.is_flat():
					# This shouldn't happen.
					raise Exception("Error when checking if flat tet is admissible")
				elif new_0.is_flat():
					if not (new_1.tilt(V3) + new_2.tilt(V2)).evaluate() > 0:
						# Then the edge is not concave, so the tet is not admissible.
						return False
				else:
					if not (new_1.tilt(V2) + new_0.tilt(V3)).evaluate() > 0:
						# Then the edge is not concave, so the tet is not admissible.
						return False
				# If we haven't returned False, then that flat tet is admissible. 
				return True
		# If we hit this, then this is a flat tet we weren't even able to do
		# a special_four_to_four move on, so it's definitely not admissible.
		return False

	"""
	PACHNER MOVES.
	"""

	"""
	This function checks if a 2-3 move through two_subsimplex of tet is possible.
	"""
	def check_two_to_three(self,two_subsimplex,tet):
		if tet.Neighbor[two_subsimplex] is None:
			return 0
		if tet.Neighbor[two_subsimplex] == tet and not tet.face_glued_to_self(two_subsimplex):
			return 0
		if self.is_geometric:
			for one_subsimplex in OneSubsimplices:
				if is_subset(one_subsimplex,two_subsimplex):
					z = tet.edge_params[one_subsimplex]
					w = tet.Neighbor[two_subsimplex].edge_params[tet.Gluing[two_subsimplex].image(one_subsimplex)]
					if (z*w).imag.evaluate() < 0:
						return 0
		for sym in tet.Symmetries:
			if sym.image(two_subsimplex) != two_subsimplex:
				return 0
		for sym in tet.Neighbor[two_subsimplex].Symmetries:
			if sym.image(tet.Gluing[two_subsimplex].image(two_subsimplex)) != tet.Gluing[two_subsimplex].image(two_subsimplex):
				return 0
		return 1

	

	"""
	The 2-3 move. Note that if build is set to 0, we do not re-build edge classes, vertex classes, and horotriangle data. We use
	this move with build = 0 in the 3-6 move, but normally the default, build = 1, should be used instead.
	"""
	def two_to_three(self, two_subsimplex, tet, build = 1):
		if tet.Neighbor[two_subsimplex] is None:
			return 0
		if self.is_geometric:
			One = self.complex_arithmetic_type.One()
		if tet.face_glued_to_self(two_subsimplex):
			for one_subsimplex in OneSubsimplices:
				if tet.Gluing[two_subsimplex].image(one_subsimplex) == one_subsimplex:
					if is_subset(one_subsimplex,two_subsimplex):
						a = Arrow(one_subsimplex,two_subsimplex,tet)
						break
		else:
			a = Arrow(PickAnEdge[two_subsimplex], two_subsimplex, tet)
			b = a.glued()
		if self.is_geometric:
			z = a.Tetrahedron.edge_params[a.south_head()]
			if tet.face_glued_to_self(two_subsimplex):
				w = z
			else:
				w = b.Tetrahedron.edge_params[b.north_tail()]
		#Now we make the new tets. We consider each case separately. That means some redundant
		#lines are written, but I think it's easier to understand this way. 
		if tet.face_glued_to_self(two_subsimplex) and tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			# Edge params.
			if self.is_geometric:
				new.Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			# Symmetries.
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new.add_sym(new.copy().reverse())
			# Face gluings.
			new.glue(new.copy().reverse())
			new.copy().opposite().true_glue_as(a.copy().reverse())
			# Link the vertices of the new tetrahedron to the pre-existing vertex classes.
			new.Tetrahedron.Class[new.north_vertex()] = a.Tetrahedron.Class[a.west_vertex()]
			new.Tetrahedron.Class[new.south_vertex()] = a.Tetrahedron.Class[a.west_vertex()]
			new.Tetrahedron.Class[new.east_vertex()] = a.Tetrahedron.Class[a.south_vertex()]
			new.Tetrahedron.Class[new.west_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
			# Assign edge labels to the new tetrahedron.
			new.Tetrahedron.edge_labels = {
				new.equator(): tet.edge_labels[a.axis()], 
				new.axis(): 3,
				new.north_head(): tet.edge_labels[a.south_tail()], 
				new.north_tail(): tet.edge_labels[a.north_tail()],
				new.south_head(): tet.edge_labels[a.north_tail()],
				new.south_tail(): tet.edge_labels[a.south_tail()]
				}
		elif tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			# Edge params.
			if self.is_geometric:
				new.Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			# Symmetries.
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			# Face gluings.
			new.glue(new.copy())
			new.copy().opposite().true_glue_as(a.copy().reverse())
			new.copy().opposite().reverse().true_glue_as(b)
			# Link the vertices of the new tetrahedron the pre-existing vertex classes.
			new.Tetrahedron.Class[new.north_vertex()] = a.Tetrahedron.Class[a.west_vertex()]
			new.Tetrahedron.Class[new.south_vertex()] = b.Tetrahedron.Class[b.east_vertex()]
			new.Tetrahedron.Class[new.east_vertex()] = a.Tetrahedron.Class[a.south_vertex()]
			new.Tetrahedron.Class[new.west_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
			# Update the edge labels. Changing arrow positions first, just so I can copy paste
			# from another time I did this.
			new.opposite().reverse()
			a.reverse()
			new.Tetrahedron.edge_labels = {
				new.equator(): 3,
				new.axis(): tet.edge_labels[a.axis()],
				new.north_head(): b.Tetrahedron.edge_labels[b.north_head()],
				new.north_tail(): tet.edge_labels[a.south_head()],
				new.south_head(): b.Tetrahedron.edge_labels[b.south_head()],
				new.south_tail(): tet.edge_labels[a.north_head()] 
				}
			# Update the canonize info.
			if tet.canonize_info is not None:
				new.Tetrahedron.canonize_info = CanonizeInfo()
				new.Tetrahedron.canonize_info.part_of_coned_cell = True
				new.Tetrahedron.canonize_info.is_flat = False
				new.Tetrahedron.canonize_info.face_status = {
					new.north_face(): 2,
					new.south_face(): 2,
					new.east_face(): b.Tetrahedron.true_face_status(b.east_face()),
					new.west_face(): a.Tetrahedron.true_face_status(a.east_face())
					}
		elif tet.face_glued_to_self(two_subsimplex):
			new = self.new_arrows(2)
			# Edge params.
			if self.is_geometric:
				new[0].Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
				new[1].Tetrahedron.fill_edge_params(z/(One-w))
			# Symmetries.
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[0].add_sym(new[0].copy().reverse())
			# Face gluings.
			new[0].glue(new[1])
			new[1].glue(new[1].copy().reverse())
			new[0].copy().opposite().true_glue_as(a.copy().reverse())
			new[1].copy().opposite().true_glue_as(a.copy().reverse().rotate(-1))
			new[1].copy().opposite().reverse().true_glue_as(a.copy().reverse().rotate(1))
			# Link the vertices of the new tetrahedra to the pre-existing vertex classes.
			for c in new:
				c.Tetrahedron.Class[c.north_vertex()] = a.Tetrahedron.Class[a.west_vertex()]
				c.Tetrahedron.Class[c.south_vertex()] = a.Tetrahedron.Class[a.west_vertex()]
				c.Tetrahedron.Class[c.east_vertex()] = a.Tetrahedron.Class[a.south_vertex()]
				c.Tetrahedron.Class[c.west_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
				a.opposite().rotate(-1)
			# Update the edge labels. Re-arranging arrows so I can copy paste this from another place I did it.
			a.reverse().rotate(1)
			new[0].opposite()
			new[1].opposite().reverse()
			new[0].Tetrahedron.edge_labels = {
				new[0].equator(): 1,
				new[0].axis(): tet.edge_labels[a.north_tail()],
				new[0].north_head(): tet.edge_labels[a.equator()],
				new[0].north_tail(): tet.edge_labels[a.north_head()],
				new[0].south_head(): tet.edge_labels[a.north_head()],
				new[0].south_tail(): tet.edge_labels[a.equator()]
				}
			new[1].Tetrahedron.edge_labels = {
				new[1].equator(): 1,
				new[1].axis(): tet.edge_labels[a.south_tail()],
				new[1].north_head(): tet.edge_labels[a.north_head()],
				new[1].north_tail(): tet.edge_labels[a.equator()],
				new[1].south_head(): tet.edge_labels[a.south_head()],
				new[1].south_tail(): tet.edge_labels[a.south_head()]
				}
		else:
			#In this case there is no face rotation nor is the face glued to itself.
			new = self.new_arrows(3)
			# Edge params.
			if self.is_geometric:
				new[0].Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
				new[1].Tetrahedron.fill_edge_params(z/(One-w))
				new[2].Tetrahedron.fill_edge_params(w/(One-z))
			# Face gluings and symmetries.
			for i in range(3):
				new[i].glue(new[(i+1)%3])
			a.reverse()
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
				c.copy().opposite().true_glue_as(a)
				c.copy().opposite().reverse().true_glue_as(b)
				a.rotate(-1)
				b.rotate(1)
			# Link the vertices of the new tetrahedra to the pre-existing vertex classes.
			for c in new:
				c.Tetrahedron.Class[c.north_vertex()] = a.Tetrahedron.Class[a.east_vertex()]
				c.Tetrahedron.Class[c.south_vertex()] = b.Tetrahedron.Class[b.east_vertex()]
				c.Tetrahedron.Class[c.east_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
				c.Tetrahedron.Class[c.west_vertex()] = a.Tetrahedron.Class[a.south_vertex()]
				a.rotate(-1)
				b.rotate(1)
			# Update the edge labels.
			for c in new:
				c.Tetrahedron.edge_labels = {
					c.equator(): tet.edge_labels[a.axis()],
					c.axis(): 1,
					c.north_head(): tet.edge_labels[a.north_head()],
					c.north_tail(): tet.edge_labels[a.south_head()],
					c.south_head(): b.Tetrahedron.edge_labels[b.south_head()],
					c.south_tail(): b.Tetrahedron.edge_labels[b.north_head()]
					}
				a.rotate(-1)
				b.rotate(1)
			# Update the canonize info.
			if tet.canonize_info is not None:	
				for c in new:
					c.Tetrahedron.canonize_info = CanonizeInfo()
					c.Tetrahedron.canonize_info.part_of_coned_cell = True
					c.Tetrahedron.canonize_info.is_flat = False
					c.Tetrahedron.canonize_info.face_status = {
						c.north_face(): a.Tetrahedron.canonize_info.face_status[a.east_face()],
						c.south_face(): b.Tetrahedron.canonize_info.face_status[b.east_face()],
						c.east_face(): 2,
						c.west_face(): 2
						}
					a.rotate(-1)
					b.rotate(1)
		self.Tetrahedra.remove(tet)
		if not tet.face_glued_to_self(two_subsimplex):
			self.Tetrahedra.remove(b.Tetrahedron)
		if build:
			self.clear_and_rebuild()
		return 1

	"""
	3-2 move. Returns 0 if a 3-2 move is not possible, otherwise does the move on self and returns 1.
	Currently it does not reverse a flat 2-3 move. Optional "build" parameter determines if we build
	all edge classes, horotriangles, etc. after doing the move. Normally we do want to build, but we 
	don't want to build when we use the 3-2 move internally in the 6-3 move.
	"""
	def three_to_two(self,edge,build = 1):
		One = ComplexSquareRootCombination.One()
		a = edge.get_arrow()
		face_glued_to_self = False
		face_rotation = False
		if edge.valence() == 1 and edge.LocusOrder == 3:
			face_rotation = True
			if len(a.Tetrahedron.Symmetries) > 2:
				return 0
			for sym in a.Tetrahedron.Symmetries:
				if sym.tuple() != (0,1,2,3):
					if sym.image(a.Edge) == a.Edge:
						face_glued_to_self = True
					else:
						return 0
		elif edge.valence() == 3 and edge.LocusOrder == 1:
			for corner in edge.Corners:
				for sym in a.Tetrahedron.Symmetries:
					if sym.image(a.Edge) != a.Edge:
						return 0
					if sym.tuple() != (0,1,2,3):
						face_glued_to_self = True
				if face_glued_to_self:
					break
				a.true_next()
			#Now, if face_glued_to_self, a is in the tet with the symmetry.
			if face_glued_to_self:
				if a.glued().Tetrahedron is None:
					a.reverse()
				check = a.copy().next()
				if len(check.Tetrahedron.Symmetries) > 1:
					return 0
				if check.next().next().reverse() != a:
					return 0
		else:
			return 0
		#If we haven't returned 0 by now, then a 3-2 move should be possible.
		if face_glued_to_self and face_rotation:
			b = self.new_arrow()
			b_copy = b.copy()
			b.glue(b_copy.reverse())
			b.reverse()
			b.add_sym(b_copy)
			b.add_sym(b_copy.rotate(1))
			b.add_sym(b_copy.rotate(1))
			v0 = a.Tetrahedron.edge_params[a.Edge]
			v1 = a.Tetrahedron.edge_params[a.Edge]
			b.Tetrahedron.fill_edge_params(One - (v0 - One)/(v0*v1 - One))
			a.opposite()
			b.true_glue_as(a)
			# Link the vertices of the new tetrahedron to the pre-existing vertex classes.
			b.Tetrahedron.Class[b.north_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
			b.Tetrahedron.Class[b.south_vertex()] = a.Tetrahedron.Class[a.south_vertex()]
			b.Tetrahedron.Class[b.east_vertex()] = a.Tetrahedron.Class[a.east_vertex()]
			b.Tetrahedron.Class[b.west_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
		elif face_rotation:
			b = self.new_arrow()
			c = self.new_arrow()
			b.glue(c)
			b.reverse()
			b_copy = b.copy()
			b.add_sym(b_copy)
			b.add_sym(b_copy.rotate(1))
			b.add_sym(b_copy.rotate(1))
			c_copy = c.copy()
			c.add_sym(c_copy)
			c.add_sym(c_copy.rotate(1))
			c.add_sym(c_copy.rotate(1))
			v0 = a.Tetrahedron.edge_params[a.Edge]
			v1 = a.Tetrahedron.edge_params[a.Edge]
			b.Tetrahedron.fill_edge_params(One - (v0 - One)/(v0*v1 - One))
			c.Tetrahedron.fill_edge_params(One - v1*(v0 - One)/(One - v1))
			a.opposite()
			b.true_glue_as(a)
			a.reverse()
			c.true_glue_as(a)
			# Link the vertices of the new tetrahedra to the pre-existing vertex classes.
			for d in [c,b]:
				d.Tetrahedron.Class[d.north_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
				d.Tetrahedron.Class[d.south_vertex()] = a.Tetrahedron.Class[a.south_vertex()]
				d.Tetrahedron.Class[d.east_vertex()] = a.Tetrahedron.Class[a.east_vertex()]
				d.Tetrahedron.Class[d.west_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
				a.reverse()
		elif face_glued_to_self:
			b = self.new_arrow()
			b.add_sym(b.copy())
			b.glue(b.copy().reverse())
			v0 = a.Tetrahedron.edge_params[a.Edge]
			v1 = a.glued().Tetrahedron.edge_params[a.glued().Edge]
			b.Tetrahedron.fill_edge_params(One - (v0 - One)/(v0*v1 - One))
			b.reverse()
			a.opposite()
			b.true_glue_as(a)
			c = a.copy().opposite().true_next()
			b.rotate(-1)
			b.glue_as(c.copy().opposite())
			b.rotate(-1)
			b.glue_as(c.copy().opposite().reverse())
			b.rotate(-1)
			# Link the vertices of the new tetrahedron to the pre-existing vertex classes.
			b.Tetrahedron.Class[b.north_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
			b.Tetrahedron.Class[b.south_vertex()] = a.Tetrahedron.Class[a.south_vertex()]
			b.Tetrahedron.Class[b.east_vertex()] = a.Tetrahedron.Class[a.east_vertex()]
			b.Tetrahedron.Class[b.west_vertex()] = c.Tetrahedron.Class[c.east_vertex()]
		else:
			b = self.new_arrow()
			c = self.new_arrow()
			b.add_sym(b.copy())
			c.add_sym(c.copy())
			v0 = a.Tetrahedron.edge_params[a.Edge]
			v1 = a.glued().Tetrahedron.edge_params[a.glued().Edge]
			b.Tetrahedron.fill_edge_params(One - (v0 - One)/(v0*v1 - One))
			c.Tetrahedron.fill_edge_params(One - v1*(v0 - One)/(One - v1))
			b.glue(c)
			b.reverse()
			# Face gluings.
			for i in range(3):
				a.opposite()
				b.glue_as(a)
				a.reverse()
				c.glue_as(a)
				b.rotate(-1)
				c.rotate(1)
				a.reverse().opposite().next()
			# Link the vertices of the new tetrahedra to pre-existing vertex classes.
			d = a.copy().next()
			a.opposite()
			for e in [b,c]:
				e.Tetrahedron.Class[e.north_vertex()] = a.Tetrahedron.Class[a.north_vertex()]
				e.Tetrahedron.Class[e.south_vertex()] = a.Tetrahedron.Class[a.south_vertex()]
				e.Tetrahedron.Class[e.east_vertex()] = a.Tetrahedron.Class[a.east_vertex()]
				e.Tetrahedron.Class[e.west_vertex()] = d.Tetrahedron.Class[d.east_vertex()]
				a.reverse()
				d.next().reverse()
		for corner in edge.Corners:
			if corner.Tetrahedron in self.Tetrahedra:
				self.Tetrahedra.remove(corner.Tetrahedron)
		if build:
			self.clear_and_rebuild()
		return 1

		
	"""
	3-6 move. Assume that tet has a single non-trivial symmetry, taking two_subsimplex to another face.
	Then the tetrahedron attached to tet along two_subsimplex (assume this is not tet) is also attached to tet 
	along that other face, making three tetrahedra in total. We can divide those three tetrahedra into six 
	which are invariant w.r.t. the symmetry. Since the symmetry identifies some of these, we only need to 
	take four of the six.

	First, we check if the 3-6 move is possible. That means we check that the symmetry is there, and
	we check certain geometric conditions are satisfied.

	Then we actually accomplish the 3-6 move by doing 3 2-3 moves. Each one on its own is illegal for
	the orbifold, but after all are done it's okay. Depending on the geometry, there are two ways to do
	the 3-6 move.

	We allow the creation of flat tetrahedra.
	"""
	def three_to_six(self,two_subsimplex,tet):
		One = ComplexSquareRootCombination.One()
		if tet.Neighbor[two_subsimplex] is None:
			return 0
		if len(tet.Symmetries) != 2:
			return 0
		sym = tet.nontrivial_sym()
		#now sym is the non-trivial symmetry.
		for one_subsimplex in OneSubsimplices:
			if is_subset(one_subsimplex,two_subsimplex) and sym.image(one_subsimplex) == one_subsimplex:
				break
		#now one_subsimplex is the edge of tet in two_subsimplex fixed by sym.
		voisin = tet.Neighbor[two_subsimplex]
		if voisin == tet or len(voisin.Symmetries) > 1:
			return 0
		a = Arrow(one_subsimplex,two_subsimplex,tet)
		b = a.glued()
		z = tet.edge_params[a.Edge]
		w = voisin.edge_params[b.Edge]
		if ((z*w*w).imag).evaluate() < 0:
			return 0
		flat_u0_u1_v1_w1 = False
		flat_u0_u1_v0_w1 = False
		z = tet.edge_params[a.south_head()]
		w = voisin.edge_params[b.south_tail()]
		if ((z*w).imag).evaluate() < 0:
			return 0
		if (z*w).imag == SquareRootCombination.Zero():
			flat_u0_u1_v1_w1 = True
		z = tet.edge_params[a.north_head()]
		w = voisin.edge_params[b.north_tail()]
		if ((z*w).imag).evaluate() < 0:
			return 0
		if (z*w).imag == SquareRootCombination.Zero():
			flat_u0_u1_v0_w1 = True
		if flat_u0_u1_v0_w1 and flat_u0_u1_v1_w1:
			#Maybe this is actually a valid case, and it could be a valid case for
			#the 2-3 move situation too. But for now, let's not allow it.
			return 0
		else:
			# Get some complex numbers which we use to determine whether [u_0,w_1] is beneath, above, or
			# intersecting [u_1,w_0].
			x = tet.edge_params[a.north_tail()]
			y = voisin.edge_params[b.Edge]
			z = voisin.edge_params[b.south_tail()]
			w0 = One - One/(y - x*y + x)
			complex_abs_w0 = ComplexSquareRootCombination(abs(w0),SquareRootCombination.Zero())
			w1 = One - One/(x*z)
			complex_abs_w1 = ComplexSquareRootCombination(abs(w1),SquareRootCombination.Zero())
			# Create a new arrow for the new tetrahedron which is the image under the symmetry of b.Tetrahedron.
			c = self.new_arrow()
			c.glue(a)
			# Link the vertices of c.Tetrahedron to the correct vertex classes.
			c.Tetrahedron.Class[c.north_vertex()] = b.Tetrahedron.Class[b.north_vertex()]
			c.Tetrahedron.Class[c.south_vertex()] = b.Tetrahedron.Class[b.south_vertex()]
			c.Tetrahedron.Class[c.east_vertex()] = a.Tetrahedron.Class[a.west_vertex()]
			c.Tetrahedron.Class[c.west_vertex()] = b.Tetrahedron.Class[b.east_vertex()]
			# Give c.Tetrahedron the geometry making it a copy of b.Tetrahedron under the symmetry.
			c.Tetrahedron.fill_edge_params(b.Tetrahedron.edge_params[b.Edge])
			if ((w0/complex_abs_w0).real.evaluate() < (w1/complex_abs_w1).real.evaluate() or
				(w0/complex_abs_w0).real == (w1/complex_abs_w1).real):
				#This is the case [u_0,w_1] "beneath" or intersecting [u_1,w_0].
				self.two_to_three(two_subsimplex,tet,0)
				c.next()
				new_c = c.copy().reverse()
				c.opposite().next().reverse()
				self.two_to_three(new_c.Face,new_c.Tetrahedron,0)
				new_c = c.copy()
				c.reverse().next()
				self.two_to_three(new_c.Face,new_c.Tetrahedron,0)
				d = c.copy()
				c.reverse().next()
				e = c.copy().opposite().reverse()
				c.add_sym(e)
				e.reverse().next()
				c.next()
				f = c.copy()
				c.opposite().reverse().next()
				c.add_sym(c.copy().opposite().reverse())
				c.rotate(-1).next()
				d.glue(e.glued())
				e.rotate(1)
				d.rotate(1)
				if e.glued().Tetrahedron != None:
					d.glue_as(e)
				e.Tetrahedron.erase()
				c.reverse().rotate(-1)
				c.glue(f.glued())
				f.Tetrahedron.erase()
				self.Tetrahedra.remove(e.Tetrahedron)
				self.Tetrahedra.remove(f.Tetrahedron)
			elif (w0/complex_abs_w0).real.evaluate() > (w1/complex_abs_w1).real.evaluate():
				#This is the case [u_0,w_1] "above" [u_1,w_0].
				self.two_to_three(c.Face,c.Tetrahedron,0)
				b.reverse()
				new_b = b.copy()
				b.next().opposite().next()
				self.two_to_three(new_b.Face,new_b.Tetrahedron,0)
				new_b = b.copy().next()
				self.two_to_three(new_b.Face,new_b.Tetrahedron,0)
				#Now add symmetries.
				d = b.copy()
				b.next()
				b.add_sym(b.copy().opposite())
				e = b.copy().opposite().reverse().next().reverse()
				b.next()
				f = b.copy()
				b.opposite().next()
				b.add_sym(b.copy().opposite())
				b.rotate(1).next()
				#Now we adjust face gluings and remove two tets.
				d.rotate(1)
				e.rotate(1)
				if d.glued().Tetrahedron != None:
					e.glue_as(d)
				d.reverse()
				e.reverse()
				e.glue(d.glued())
				d.Tetrahedron.erase()
				b.rotate(-1)
				f.reverse().rotate(1).reverse()
				f.glue(b.glued())
				b.Tetrahedron.erase()
				self.Tetrahedra.remove(d.Tetrahedron)
				self.Tetrahedra.remove(b.Tetrahedron)
		for tet in self.Tetrahedra:
			tet.fix_glued_to_self()
		self.clear_and_rebuild()
		return 1

	# Reverse of 3-6. This hasn't really been used yet, it might have bugs. It seems like it's
	# not actually necessary for canonize.
	def six_to_three(self,edge):
		#Try to do a 6-3 move, viewing edge as [w0,w1]. This gets complicated and requires a lot
		#of checking before you know you can do the move. See my write-up.
		if edge.valence() != 4 or edge.LocusOrder != 1:
			return 0
		for corner in edge.Corners:
			if len(corner.Tetrahedron.Symmetries) == 1:
				middle = Arrow(corner.Subsimplex,RightFace[corner.Subsimplex],corner.Tetrahedron)
				break
		else:
			return 0
		a = middle.glued()
		b = middle.reverse().glued()
		if a.Tetrahedron is b.Tetrahedron:
			return 0
		for c in [a,b]:
			if len(c.Tetrahedron.Symmetries) != 2:
				return 0
			for sym in c.Tetrahedron.Symmetries:
				if sym.image(c.Edge) != c.Edge:
					return 0
		#This is the end of part 1 of our check.
		setup_success = False
		for c in [b,a]:
			#This is not a typo. I do want to start with b.
			d = c.copy().reverse().rotate(1)
			e = d.copy()
			if e.glued() is None:
				e.opposite().reverse()
			e.next()
			if e.Tetrahedron is middle.Tetrahedron or len(e.Tetrahedron.Symmetries) != 1:
				continue
			if e.glued() == middle.copy().opposite().rotate(-1):
				#Then a 6-3 move is possible and we're in case 1.
				case = 1
				setup_success = True
				break
			e.opposite()
			if e.glued() == middle.copy().opposite().reverse():
				case = 2
				setup_success = True
				break
			#If we got here with c = b, then it doesn't work with b as bottom. Now we want to
			#check if a could be bottom. The check is exactly the same computation as long as we
			#reverse middle before. So we do that, then re-enter the loop.
			middle.reverse()
		#Now we fix some labels and adjust face gluings of d.Tetrahedron if the wrong face is glued to None.
		if setup_success:
			bottom = c
			if c is a:
				top = b
			else:
				top = a
			fourth = e
			f = d.copy().opposite().reverse()
			if case == 1 and d.glued() is None:
				d.glue(f.glued())
				f.Tetrahedron.Neighbor[f.Face] = None
				f.Tetrahedron.Gluing[f.Face] = None
			if case == 2 and d.glued() is not None:
				f.glue(d.glued())
				d.Tetrahedron.Neighbor[d.Face] = None
				d.Tetrahedron.Gluing[d.Face] = None
		else:
			return 0
		#We are finally done with all the checking. If we got to this point, a 6-3 move is possible.
		#Now we just reverse the 2-3 moves done in the 3-6 move, which depend on the case. But first
		#we need to add two tets, which are redundant by the symmetries but needed for this move.
		new_middle_tet = Tetrahedron()
		new_middle_tet.fill_edge_params(middle.Tetrahedron.edge_params[E01])
		new_middle_tet.Symmetries.append(Perm4((0,1,2,3)))
		new_middle = Arrow(middle.Edge,middle.Face,new_middle_tet)
		top.glue(new_middle)
		bottom.glue(new_middle.copy().reverse())
		new_fourth_tet = Tetrahedron()
		new_fourth_tet.fill_edge_params(fourth.Tetrahedron.edge_params[E01])
		new_fourth_tet.Symmetries.append(Perm4((0,1,2,3)))
		new_fourth = Arrow(fourth.Edge,fourth.Face,new_fourth_tet)
		self.add_tet(new_middle_tet)
		self.add_tet(new_fourth_tet)
		#Now we carefully glue and re-glue faces so the union of all six tets is the particular
		#octahedron as in the main picture.
		if case == 1:
			#face [u1,v0,w0] of middle.Tet is now glued to new_fourth_tet.
			new_fourth.glue(middle.copy().opposite().rotate(-1))
			#face [u0,v1,w1] of fourth.Tet is now glued to new_middle_tet.
			fourth.glue(new_middle.copy().opposite().rotate(-1))
			#face [u0,u1,w0] of new_fourth_tet is glued to bottom.tet
			new_fourth.copy().reverse().glue(bottom.copy().rotate(1).reverse())
			#The two remaining faces of new_fourth_tet, [u0,u1,v0] and [u0,v0,w0],
			#remain glued to None. They're external faces in the ocahedron which
			#should be glued to None.
		if case == 2:
			#face [u1,v1,w0] of middle.Tet is now glued to new_fourth_tet.
			new_fourth.glue(middle.copy().opposite().reverse())
			#face [u0,v0,w1] of fourth.Tet is now glued to new_middle_tet.
			fourth.glue(new_middle.copy().opposite().reverse())
			#face [u0,u1,w0] of new_fourth_tet is glued to bottom.Tet.
			new_fourth.copy().rotate(1).glue(bottom.copy().opposite().reverse())
			#The two remaining faces of new_fourth_tet will remain glued to None.
		#We're going to do some 3-2 moves to reverse the 2-3 moves which are in the 3-6 move.
		#before we do that, we should remove the nontrivial symmetries from any of our six tets
		#which have them. This is because the symmetries will confuse the three_to_two function.
		top.Tetrahedron.Symmetries = [Perm4((0,1,2,3))]
		bottom.Tetrahedron.Symmetries = [Perm4((0,1,2,3))]
		#Now we can do the 3-2 moves. The main argument of the three_to_two function is an
		#edge class. At this point, we've messed with the triangulation so much that the
		#edges in self.Edges are meaningless. So we should make a new, temporary, edge for
		#each 3-2 move we do. In fact we use the same Edge object for each 3-2 move, but clear
		#its Corners after each use.
		temp_edge = Edge()
		temp_edge.LocusOrder = 1 
		if case == 1:
			#Define the arrow a to help us save an arrow in one of the newly created tetrahedra,
			#i.e. one of the three tets resulting after all 3-2 moves. We will then use a to
			#manipulate them as needed.
			a = fourth.copy() 
			arrows = [new_fourth,top,fourth]
			for i in range(3):
				if i == 2:
					a.next().opposite().reverse().next()
					#Now a is in [u0,v0,v1,w0] pointing from u0 to w0. 
				for j in range(3):
					temp_edge.Corners.append(Corner(arrows[i].Tetrahedron,arrows[i].Edge))
					arrows[i].next()
				self.three_to_two(temp_edge,0)
				#The 0 parameter is used so we don't try to build edge classes, horotriangles, etc. 
				#after applying the 3-2 move.
				temp_edge.Corners = []
			a.reverse()
			a.glued().add_sym(a.glued().copy().reverse())
			a.Tetrahedron.detach(a.Face)
			self.Tetrahedra.remove(a.Tetrahedron)
		if case == 2:
			new_fourth.opposite().rotate(-1)
			#Define an arrow a which will be used in the same sense here as in case 1.
			a = new_fourth.copy()
			arrows = [bottom.reverse().rotate(-1),top,new_fourth]
			for i in range(3):
				if i == 2:
					a.next().opposite().next()
					#now a is in [u1,v0,v1,w1], pointing from u1 to w1.
				for j in range(3):
					temp_edge.Corners.append(Corner(arrows[i].Tetrahedron,arrows[i].Edge))
					arrows[i].next()
				self.three_to_two(temp_edge,0)
				temp_edge.Corners = []
			a.reverse().next()
			a.add_sym(a.copy().reverse())
			a.next().reverse()
			a.Tetrahedron.detach(a.Face)
			self.Tetrahedra.remove(a.Tetrahedron)
		self.clear_and_rebuild()
		return 1


	"""
	4-4 move around a valence 4 edge. In the case that there are 4 distinct tetrahedra, none of them
	having symmetries, this is just the usual manifold 4-4 move. 

	The union of the 4 tetrahedra is an octahedron. We might consider an orbifold version of the 4-4
	move for each group of symmetries of an octahedron. For my purposes, I'm only considering one such
	case. The case where you have an order 2 symmetry, pi rotation through an axis connecting opposite
	edges of the octahedron. See the write-up for a picture and explanation.

	UPDATE: To reverse the 4-4 with symmetry, use the special_four_to_four function.
	
	So the symmetry can only be pi rotation swapping x1 with x4 and x2 with x3,
	where the labels are from my write-up. 
	"""	
	def four_to_four(self,edge):
		if edge.valence() != 4:
			return 0
		for corner in edge.Corners:
			if len(corner.Tetrahedron.Symmetries) == 1:
				a = Arrow(corner.Subsimplex,RightFace[corner.Subsimplex],corner.Tetrahedron)
				break
		else:
			return 0
		b = a.glued()
		c = a.reverse().glued()
		for d in [b,c]:
			for sym in d.Tetrahedron.Symmetries:
				if sym.image(d.Edge) != d.Edge:
					return 0
		if len(b.Tetrahedron.Symmetries) == 2 and len(c.Tetrahedron.Symmetries) == 2:
			symmetry = True
		elif len(b.Tetrahedron.Symmetries) == 1 and len(c.Tetrahedron.Symmetries) == 1:
			symmetry = False
		else:
			return 0
		za = a.Tetrahedron.edge_params[a.north_head()]
		zb = b.Tetrahedron.edge_params[b.Edge]
		zc = c.Tetrahedron.edge_params[c.north_tail()]
		One = ComplexSquareRootCombination.One()
		wa = za
		wb = za*zb/(zb + za - One)
		wc = za*zc
		if not symmetry:
			old = [c.next() for i in range(4)]
			new = self.new_arrows(4)
			new[0].Tetrahedron.fill_edge_params(wc*(wb - One)/((wc - One)*wb))
			new[1].Tetrahedron.fill_edge_params(wc)
			new[2].Tetrahedron.fill_edge_params((wa - wc)/(One - wc))
			new[3].Tetrahedron.fill_edge_params((wa - wc)*(wb - One)/((wa - One)*(wb - wc)))
			new[0].reverse().rotate(1)
			new[1].opposite().reverse()
			new[2].rotate(1).reverse()
			#new[3] is already in the right place.
			# Now the old and new arrows are in position as in the diagram. We now
			# link the vertices of the new tetrahedra to the correct vertex classes.
			# new[0]
			new[0].Tetrahedron.Class[new[0].north_vertex()] = old[1].Tetrahedron.Class[old[1].east_vertex()]
			new[0].Tetrahedron.Class[new[0].south_vertex()] = old[0].Tetrahedron.Class[old[0].west_vertex()]
			new[0].Tetrahedron.Class[new[0].east_vertex()] = old[1].Tetrahedron.Class[old[1].west_vertex()]
			new[0].Tetrahedron.Class[new[0].west_vertex()] = old[1].Tetrahedron.Class[old[1].north_vertex()]
			# new[1]
			new[1].Tetrahedron.Class[new[1].north_vertex()] = old[1].Tetrahedron.Class[old[1].east_vertex()]
			new[1].Tetrahedron.Class[new[1].south_vertex()] = old[0].Tetrahedron.Class[old[0].west_vertex()]
			new[1].Tetrahedron.Class[new[1].east_vertex()] = old[1].Tetrahedron.Class[old[1].south_vertex()]
			new[1].Tetrahedron.Class[new[1].west_vertex()] = old[1].Tetrahedron.Class[old[1].west_vertex()]
			# new[2]
			new[2].Tetrahedron.Class[new[2].north_vertex()] = old[1].Tetrahedron.Class[old[1].east_vertex()]
			new[2].Tetrahedron.Class[new[2].south_vertex()] = old[0].Tetrahedron.Class[old[0].west_vertex()]
			new[2].Tetrahedron.Class[new[2].east_vertex()] = old[2].Tetrahedron.Class[old[2].east_vertex()]
			new[2].Tetrahedron.Class[new[2].west_vertex()] = old[1].Tetrahedron.Class[old[1].south_vertex()]
			# new[3]
			new[3].Tetrahedron.Class[new[3].north_vertex()] = old[1].Tetrahedron.Class[old[1].east_vertex()]
			new[3].Tetrahedron.Class[new[3].south_vertex()] = old[0].Tetrahedron.Class[old[0].west_vertex()]
			new[3].Tetrahedron.Class[new[3].east_vertex()] = old[1].Tetrahedron.Class[old[1].north_vertex()]
			new[3].Tetrahedron.Class[new[3].west_vertex()] = old[2].Tetrahedron.Class[old[2].east_vertex()]
			# Now we do the gluings.
			for i in range(4):
				new[i].glue(new[(i+1)%4])
			#gluing faces of new[0]
			old[0].opposite()
			new[0].reverse().rotate(-1)
			new[0].glue_as(old[0])
			old[0].opposite()
			new[0].rotate(-1)
			old[1].opposite()
			new[0].glue_as(old[1])
			old[1].opposite()
			#now glue faces of new[1]
			new[1].rotate(1)
			old[0].opposite().reverse()
			new[1].glue_as(old[0])
			new[1].rotate(1)
			old[1].opposite().reverse()
			new[1].glue_as(old[1])
			#now glue faces of new[2]
			new[2].reverse().rotate(1)
			old[2].opposite().reverse()
			new[2].glue_as(old[2])
			old[2].opposite().reverse()
			new[2].rotate(1)
			old[3].opposite().reverse()
			new[2].glue_as(old[3])
			old[3].opposite().reverse()
			#now faces of new[3]
			new[3].rotate(-1)
			old[2].opposite()
			new[3].glue_as(old[2])
			new[3].rotate(-1)
			old[3].opposite()
			new[3].glue_as(old[3])
			for i in range(4):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			for a in old:
				self.Tetrahedra.remove(a.Tetrahedron)
		else:
			#Because there is no tet at [x1,x2,y1,y2], have to awkwardly define old
			#to be consistent with how I chose to write this.
			old = []
			old.append(a.copy().reverse())
			old.append(b.reverse())
			old.append(a)
			old.append(c)
			new = self.new_arrows(3)
			#new[0] is [x1,x3,x4,y1]
			new[0].Tetrahedron.fill_edge_params(wc)
			#new[1] is [x1,x2,x3,y1]
			new[1].Tetrahedron.fill_edge_params(wc*(wb - One)/((wc - One)*wb))
			#new[2] is [x1,x2,x3,x4], which could be flat.
			new[2].Tetrahedron.fill_edge_params((wb - wc)/(One - wc))
			new[0].opposite().reverse()
			new[1].reverse().rotate(1)
			new[2].rotate(1).reverse()
			# Now old and new arrows are in correct positions as in the diagram.
			# Link the vertices of the new tetahedra to the correct pre-existing vertex classes.
			# new[0]
			new[0].Tetrahedron.Class[new[0].north_vertex()] = old[1].Tetrahedron.Class[old[1].east_vertex()]
			new[0].Tetrahedron.Class[new[0].south_vertex()] = old[0].Tetrahedron.Class[old[0].west_vertex()]
			new[0].Tetrahedron.Class[new[0].east_vertex()] = old[1].Tetrahedron.Class[old[1].north_vertex()]
			new[0].Tetrahedron.Class[new[0].west_vertex()] = old[2].Tetrahedron.Class[old[2].east_vertex()]
			# new[1]
			new[1].Tetrahedron.Class[new[1].north_vertex()] = old[1].Tetrahedron.Class[old[1].east_vertex()]
			new[1].Tetrahedron.Class[new[1].south_vertex()] = old[0].Tetrahedron.Class[old[0].west_vertex()]
			new[1].Tetrahedron.Class[new[1].east_vertex()] = old[1].Tetrahedron.Class[old[1].west_vertex()]
			new[1].Tetrahedron.Class[new[1].west_vertex()] = old[1].Tetrahedron.Class[old[1].north_vertex()]
			# new[2]
			new[2].Tetrahedron.Class[new[2].north_vertex()] = old[1].Tetrahedron.Class[old[1].east_vertex()]
			new[2].Tetrahedron.Class[new[2].south_vertex()] = old[0].Tetrahedron.Class[old[0].west_vertex()]
			new[2].Tetrahedron.Class[new[2].east_vertex()] = old[2].Tetrahedron.Class[old[2].east_vertex()]
			new[2].Tetrahedron.Class[new[2].west_vertex()] = old[1].Tetrahedron.Class[old[1].west_vertex()]
			# Now do the face gluings.
			for i in range(3):
				new[i].glue(new[(i+1)%3])
			#make the other face gluings for new[0]
			new[0].rotate(-1)
			old[2].opposite()
			new[0].glue_as(old[2])
			new[0].rotate(-1)
			old[3].opposite()
			new[0].true_glue_as(old[3])
			#now for new[1]
			new[1].reverse().rotate(-1)
			old[0].opposite()
			new[1].glue_as(old[0])
			new[1].rotate(-1)
			old[1].opposite()
			new[1].true_glue_as(old[1])
			#now add the symmetries.
			for i in range(3):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[2].add_sym(new[2].copy().opposite().reverse())
			for a in old:
				if a.Tetrahedron in self.Tetrahedra:
					self.Tetrahedra.remove(a.Tetrahedron)
		self.clear_and_rebuild()
		return 1

	
	"""
	This 4-4 move is the reverse of the one above in the case of a symmetry. It could work even
	if there isn't a flat tet, but I choose to have it return 0 then. The reason is that, if there
	is not a flat tet, then some other pachner move will help you make progress when finding the
	canonical decomposition. When there is a flat tet there, it's a special situation which must
	be dealt with carefully when trying to canonize. 

	The non-geometric version of this move is programmed in SimplicialOrbifold.
	All that's new here is determining geometry.

	Something that I didn't realize until recently is that there are two cases. When the symmetry axis
	is horizontal, and when it's vertical.
	"""
	def special_four_to_four(self,edge):
		# This is the case where we have half of an octahedron, with the square face glued to 
		# itself by rotation around an axis which intersects the midpoints of two sides.
		if edge.valence() != 3 or edge.LocusOrder != 1:
			return 0
		a = edge.get_arrow()
		if len(a.Tetrahedron.Symmetries) != 2:
			# Let's move the arrow into the tet with a nontrivial symmetry (if it exists).
			a.true_next()
			if len(a.Tetrahedron.Symmetries) != 2:
				a.true_next()
			if len(a.Tetrahedron.Symmetries) != 2:
				return 0
		# Now the arrow a is positioned as in the diagram (in the write-up) and we do
		# the final checks.
		if a.Tetrahedron.is_flat() == False:
			return 0
		sym = a.Tetrahedron.nontrivial_sym()
		if sym.image(a.Edge) == a.Edge:
			return 0
		b = a.copy().true_next()
		c = b.glued()
		if b.Tetrahedron is None or c.Tetrahedron is None:
			return 0
		if b.Tetrahedron is c.Tetrahedron:
			return 0
		if len(b.Tetrahedron.Symmetries) != 1 or len(c.Tetrahedron.Symmetries) != 1:
			return 0
		# If we've gotten to here, then this case of a 4-4 move can be done.
		# We must determine if the symmetry axis is horizontal or vertical.
		if comp(a.south_face()) == sym.image(a.head()):
			case = 'horizontal'
		elif comp(a.north_face()) == sym.image(a.head()):
			case = 'vertical'
		else:
			raise Exception('error in special_four_to_four')
		new = self.new_arrows(3)
		# First we take care of geometry.
		One = ComplexSquareRootCombination.One()
		# u0 = shape of (y1,x4) in (y1,x1,x2,x4)
		u0 = b.Tetrahedron.edge_params[b.north_head()]
		# u1 = shape of (y1,x4) in (y1,x2,x3,x4)
		u1 = c.Tetrahedron.edge_params[c.north_tail()]
		# w0, w1, and w2 are points in the complex plane, ideal points of the octahedron.
		w0 = u0*u1
		w1 = u1
		if case == 'horizontal':
			w2 = One - One/u0 + u1
		if case == 'vertical':
			w2 = (One - u1)*(u0*u1 - One) + One
		# z0 = shape of (y1,x4) in (y1,y2,x1,x4)
		z0 = w0/w2
		# z1 = shape of (y1,y2) in (y1,y2,x2,x3)
		z1 = (One - w2)*w1/(w1 - w2)
		# z2 = shape of (y1,x4) in (y1,y2,x3,x4)
		z2 = w2
		# z3 = shape of (y1,y2) in (y1,y2,x1,x2)
		z3 = (w1 - w2)*w0/(w1*(w0 - w2))		
		if case == 'horizontal':
			new[0].Tetrahedron.fill_edge_params(z1)
			new[1].Tetrahedron.fill_edge_params(One/(One - z2))
			new[2].Tetrahedron.fill_edge_params(One - One/z0)
			for i in range(3):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			for i in range(2):
				new[i].glue(new[i+1])
			new[0].add_sym(new[0].copy().reverse())
			new[2].add_sym(new[2].copy().reverse())
			new[0].copy().opposite().glue_as(c.copy().reverse().rotate(-1))
			new[1].copy().opposite().glue_as(c.copy().reverse().rotate(1))
			new[1].copy().opposite().reverse().glue_as(b.copy().rotate(1))
			new[2].copy().opposite().glue_as(b.copy().rotate(-1))
			# Link the vertices of the new tetrahedra to existing vertex classes.
			# new[0]
			new[0].Tetrahedron.Class[new[0].north_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[0].Tetrahedron.Class[new[0].south_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[0].Tetrahedron.Class[new[0].east_vertex()] = c.Tetrahedron.Class[c.east_vertex()]
			new[0].Tetrahedron.Class[new[0].west_vertex()] = c.Tetrahedron.Class[c.south_vertex()]
			# new[1]
			new[1].Tetrahedron.Class[new[1].north_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[1].Tetrahedron.Class[new[1].south_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[1].Tetrahedron.Class[new[1].east_vertex()] = c.Tetrahedron.Class[c.north_vertex()]
			new[1].Tetrahedron.Class[new[1].west_vertex()] = c.Tetrahedron.Class[c.east_vertex()]
			# new[2]
			new[2].Tetrahedron.Class[new[2].north_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[2].Tetrahedron.Class[new[2].south_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[2].Tetrahedron.Class[new[2].east_vertex()] = b.Tetrahedron.Class[b.west_vertex()]
			new[2].Tetrahedron.Class[new[2].west_vertex()] = c.Tetrahedron.Class[c.north_vertex()]
		if case == 'vertical':
			new[0].Tetrahedron.fill_edge_params(z3)
			new[1].Tetrahedron.fill_edge_params(z1)
			new[2].Tetrahedron.fill_edge_params(One/(One - z2))
			for i in range(3):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			for i in range(2):
				new[i].glue(new[i+1])
			new[0].add_sym(new[0].copy().reverse())
			new[2].add_sym(new[2].copy().reverse())
			new[0].copy().opposite().glue_as(b.copy().rotate(1))
			new[1].copy().opposite().glue_as(c.copy().reverse().rotate(-1))
			new[1].copy().opposite().reverse().glue_as(b.copy().rotate(-1))
			new[2].copy().opposite().glue_as(c.copy().reverse().rotate(1))
			# Link the vertices of the new tetrahedra to existing vertex classes.
			# new[0]
			new[0].Tetrahedron.Class[new[0].north_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[0].Tetrahedron.Class[new[0].south_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[0].Tetrahedron.Class[new[0].east_vertex()] = c.Tetrahedron.Class[c.south_vertex()]
			new[0].Tetrahedron.Class[new[0].west_vertex()] = b.Tetrahedron.Class[b.west_vertex()]
			# new[1]
			new[1].Tetrahedron.Class[new[1].north_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[1].Tetrahedron.Class[new[1].south_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[1].Tetrahedron.Class[new[1].east_vertex()] = c.Tetrahedron.Class[c.east_vertex()]
			new[1].Tetrahedron.Class[new[1].west_vertex()] = c.Tetrahedron.Class[c.south_vertex()]
			# new[2]
			new[2].Tetrahedron.Class[new[2].north_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[2].Tetrahedron.Class[new[2].south_vertex()] = c.Tetrahedron.Class[c.west_vertex()]
			new[2].Tetrahedron.Class[new[2].east_vertex()] = c.Tetrahedron.Class[c.north_vertex()]
			new[2].Tetrahedron.Class[new[2].west_vertex()] = c.Tetrahedron.Class[c.east_vertex()]
		self.Tetrahedra.remove(a.Tetrahedron)
		self.Tetrahedra.remove(b.Tetrahedron)
		self.Tetrahedra.remove(c.Tetrahedron)
		self.clear_and_rebuild()
		return 1



	"""
	We use the following to try to get rid of flat tetrahedra. If two are glued to themselves in a
	certain way, or one is glued to itself in a certain way, then it will be "cancelled" out by this
	function. This is the orbifold version of a 2-0 move. See the write-up for more detail.
	"""
	def cancel_tetrahedra(self,edge):
		if edge.valence() == 2 and edge.LocusOrder == 1:
			#Then the adjacent tet(s) are/is definitely flat. The possible cases are:
			#1. A single tet with faces properly glued to themselves, and without the symmetry.
			#2. Two tets without the symmetry.
			#3. Two tets both with the symmetry.
			corner = edge.Corners[0]
			a = Arrow(corner.Subsimplex,LeftFace[corner.Subsimplex],corner.Tetrahedron)
			if a.copy().next() is None:
				a.reverse()
			#We don't want to try the cancellation if there are certain nontrivial symmetries.
			#They probably shouldn't be there anyway for a geometric triangulation, but we
			#check for them just in case.
			for sym in a.Tetrahedron.Symmetries:
				if sym.image(a.Edge) != a.Edge:
					return 0
			#There are two main cases. If the faces adjacent to "edge" are glued to themselves,
			#or not.
			if a.Tetrahedron.face_glued_to_self(a.Face):
				#Check it's glued to itself properly.
				perm = a.Tetrahedron.Gluing[a.Face]
				if perm.image(a.Edge) != a.Edge:
					return 0
				b = a.copy().opposite()
				if len(a.Tetrahedron.Symmetries) == 2:
					if b.copy().next() is None:
						b.reverse()
					b.next().reverse()
					b.glue(b.copy().reverse())
				else:
					#Then the only symmetry is the identity.
					c = b.copy().reverse()
					b.next().reverse()
					b.glue(c.glued())
				self.Tetrahedra.remove(a.Tetrahedron)
			else:
				#In this case there are two flat tetrahedra.
				if a.copy().next() is None:
					a.reverse()
				b = a.glued()
				#Make sure the other tet doesn't have any illegal symmetries.
				for sym in b.Tetrahedron.Symmetries:
					if sym.image(b.Edge) != b.Edge:
						return 0
				#Another sanity check.
				if len(a.Tetrahedron.Symmetries) != len(b.Tetrahedron.Symmetries):
					return 0
				#Now there are two cases: they have the nontrivial symmetry or not.
				a.opposite()
				b.opposite()
				if len(a.Tetrahedron.Symmetries) == 1:
					c = a.copy().next().reverse()
					c.glue(b.glued())
					a.reverse()
					b.reverse()
					c = a.copy().next().reverse()
					c.glue(b.glued())
				else:
					if a.copy().next() is None:
						a.reverse()
						b.reverse()
					c = a.copy().next().reverse()
					if b.copy().next() is None:
						c.glue(b.reverse().glued())
					else:
						c.glue(b.glued())
				self.Tetrahedra.remove(a.Tetrahedron)
				self.Tetrahedra.remove(b.Tetrahedron)
		elif edge.valence() == 1 and edge.LocusOrder == 2:
			#Then there's a single flat tet, with faces adjacent to the edge glued to each other. 
			#It could have the symmetry or not.
			corner = edge.Corners[0]
			a = Arrow(corner.Subsimplex,LeftFace[corner.Subsimplex],corner.Tetrahedron)
			if a.copy().next() is None:
				a.reverse()
			#We don't want to try the cancellation if there are certain nontrivial symmetries.
			#They probably shouldn't be there anyway for a geometric triangulation, but we
			#check for them just in case.
			for sym in a.Tetrahedron.Symmetries:
				if sym.image(a.Edge) != a.Edge:
					return 0
			b = a.copy().opposite()
			if len(a.Tetrahedron.Symmetries) == 2:
				if b.copy().next() is None:
					b.reverse()
				b.next().reverse()
				b.glue(b.copy().reverse())
			else:
				b.next().reverse()
				b.glue(b.copy().reverse())
				b = a.copy().opposite().reverse()
				b.next().reverse()
				b.glue(b.copy().reverse())
			self.Tetrahedra.remove(a.Tetrahedron)
		else:
			#In this case there aren't flat tets around this edge.
			return 0
		self.clear_and_rebuild()
		return 1



	def retriangulate_cube(self,two_subsimplex,tet):
		neighbor = tet.Neighbor[two_subsimplex]
		if neighbor is None or neighbor is tet:
			return 0
		if not tet.is_regular() or not neighbor.is_regular():
			return 0
		if len(tet.Symmetries) == 2:
			for sym in tet.Symmetries:
				if sym.tuple() != (0,1,2,3):
					break
			for other_face in TwoSubsimplices:
				if other_face != two_subsimplex and other_face != sym.image(two_subsimplex):
					other_neighbor = tet.Neighbor[other_face]
					if other_neighbor is not None:
						break
			if not other_neighbor.is_regular():
				return 0
			if len(neighbor.Symmetries) != 1 or len(other_neighbor.Symmetries) != 1:
				return 0
			case = 2
		elif len(tet.Symmetries) == 12:
			# Then neighbor certainly has the order 3 rotations fixing the common face.
			# We need to make sure it doesn't have any other symmetries.
			if len(neighbor.Symmetries) != 3:
				return 0
			case = 1
		elif len(tet.Symmetries) == 4:
			# Need to make sure neighbor doesn't have any non-trivial symmetries.
			# Otherwise we can't do the cube re-triangulation.
			if len(neighbor.Symmetries) != 1:
				return 0
			case = 3
		else:
			return 0
		# I'm not going to program case 2 right now. It might not be needed.
		if case == 1:
			for one_subsimplex in OneSubsimplices:
				if is_subset(one_subsimplex,two_subsimplex):
					a = Arrow(one_subsimplex,two_subsimplex,tet)
					b = a.glued()
					if b.copy().next() is not None:
						break
			#Now a and b are as in the picture in the write-up.
			new = self.new_arrows(3)
			#Tets of new[0] and new[1] are regular, new[2] is flat.
			u = SquareRootCombination([(1,Fraction(1/2))])
			v = SquareRootCombination([(Fraction(3,1),Fraction(1/2))])
			z = ComplexSquareRootCombination(u,v)
			One = ComplexSquareRootCombination.One()
			Two = ComplexSquareRootCombination(SquareRootCombination.Two(),SquareRootCombination.Zero())
			new[0].Tetrahedron.fill_edge_params(z)
			new[1].Tetrahedron.fill_edge_params(z)
			new[2].Tetrahedron.fill_edge_params(One/Two)
			new[2].reverse().rotate(1)
			new[0].glue(new[1])
			new[1].glue(new[2])
			# Now the arrows are as in the picture.
			# Let's link the vertices of the new tetrahedra to the existing vertex classes. Because of the symmetries
			# of the cube, there are only two vertex classes we have to consider.
			a_tet_vertex = a.Tetrahedron.Class[a.north_vertex()]
			b_tet_vertex = b.Tetrahedron.Class[b.east_vertex()]
			# new[0]
			new[0].Tetrahedron.Class[new[0].north_vertex()] = b_tet_vertex
			new[0].Tetrahedron.Class[new[0].south_vertex()] = b_tet_vertex
			new[0].Tetrahedron.Class[new[0].east_vertex()] = b_tet_vertex
			new[0].Tetrahedron.Class[new[0].west_vertex()] = b_tet_vertex
			# new[1]
			new[1].Tetrahedron.Class[new[1].north_vertex()] = b_tet_vertex
			new[1].Tetrahedron.Class[new[1].south_vertex()] = b_tet_vertex
			new[1].Tetrahedron.Class[new[1].east_vertex()] = a_tet_vertex
			new[1].Tetrahedron.Class[new[1].west_vertex()] = b_tet_vertex
			# new[2]
			new[2].Tetrahedron.Class[new[2].north_vertex()] = b_tet_vertex
			new[2].Tetrahedron.Class[new[2].south_vertex()] = b_tet_vertex
			new[2].Tetrahedron.Class[new[2].east_vertex()] = a_tet_vertex
			new[2].Tetrahedron.Class[new[2].west_vertex()] = a_tet_vertex
			new[0].Tetrahedron.Symmetries = [Perm4((perm)) for perm in Perm4._rawA4]
			new[1].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[1].add_sym(new[1].copy().rotate(1))
			new[1].add_sym(new[1].copy().rotate(-1))
			new[2].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[2].add_sym(new[2].copy().reverse())
			new[2].opposite().reverse()
			new[2].glue_as(b)
			self.Tetrahedra.remove(a.Tetrahedron)
			self.Tetrahedra.remove(b.Tetrahedron)
		elif case == 2:
			return 0
			#might add this later.
		elif case == 3:
			# Similar to case 1. But there are three new flat tets.
			for one_subsimplex in OneSubsimplices:
				if is_subset(one_subsimplex,two_subsimplex):
					a = Arrow(one_subsimplex,two_subsimplex,tet)
					b = a.glued()
					break
			#Now a and b are as in the picture in the write-up.
			new = self.new_arrows(5)
			#Tets of new[0] and new[1] are regular, new[2], new[3],new[4] are flat.
			u = SquareRootCombination([(1,Fraction(1/2))])
			v = SquareRootCombination([(Fraction(3,1),Fraction(1/2))])
			z = ComplexSquareRootCombination(u,v)
			One = ComplexSquareRootCombination.One()
			Two = ComplexSquareRootCombination(SquareRootCombination.Two(),SquareRootCombination.Zero())
			new[0].Tetrahedron.fill_edge_params(z)
			new[1].Tetrahedron.fill_edge_params(z)
			new[2].Tetrahedron.fill_edge_params(One/Two)
			new[3].Tetrahedron.fill_edge_params(One/Two)
			new[4].Tetrahedron.fill_edge_params(One/Two)
			new[2].reverse().rotate(1)
			new[3].reverse().rotate(1)
			new[4].reverse().rotate(1)
			new[0].glue(new[1])
			new[1].glue(new[2])
			# Now the arrows are as in the picture.
			# We now link the vertices of the new tetrahedra to the pre-existing vertex classes. As in case 1, the symmetries of
			# the cube imply its vertices have two classes. One class for the vertices of a.Tetrahedron, and the other for the
			# vertex of b.Tetrahedron opposite a.Tetrahedron. They could be the same class, depending on what's happening in the
			# rest of the triangulation, but that doesn't matter to us here.
			a_tet_vertex = a.Tetrahedron.Class[a.north_vertex()]
			b_tet_vertex = b.Tetrahedron.Class[b.east_vertex()]
			# new[0]
			new[0].Tetrahedron.Class[new[0].north_vertex()] = b_tet_vertex
			new[0].Tetrahedron.Class[new[0].south_vertex()] = b_tet_vertex
			new[0].Tetrahedron.Class[new[0].east_vertex()] = b_tet_vertex
			new[0].Tetrahedron.Class[new[0].west_vertex()] = b_tet_vertex
			# new[1]
			new[1].Tetrahedron.Class[new[1].north_vertex()] = b_tet_vertex
			new[1].Tetrahedron.Class[new[1].south_vertex()] = b_tet_vertex
			new[1].Tetrahedron.Class[new[1].east_vertex()] = a_tet_vertex
			new[1].Tetrahedron.Class[new[1].west_vertex()] = b_tet_vertex
			# new[2], new[3], new[4]
			for new in [new[2],new[3],new[4]]:
				new.Tetrahedron.Class[new.north_vertex()] = b_tet_vertex
				new.Tetrahedron.Class[new.south_vertex()] = b_tet_vertex
				new.Tetrahedron.Class[new.east_vertex()] = a_tet_vertex
				new.Tetrahedron.Class[new.west_vertex()] = a_tet_vertex
			# Symmetries and face gluings.
			new[0].Tetrahedron.Symmetries = [
				Perm4((0,1,2,3)),
				Perm4((1,0,3,2)),
				Perm4((2,3,0,1)),
				Perm4((3,2,1,0))]
			new[1].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[2].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[2].add_sym(new[2].copy().reverse())
			new[3].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[3].add_sym(new[3].copy().reverse())
			new[4].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[4].add_sym(new[4].copy().reverse())
			new[2].opposite().reverse()
			new[2].glue_as(b)
			new[1].rotate(-1)
			new[1].glue(new[3])
			new[3].opposite().reverse()
			new[3].glue_as(b.copy().rotate(1))
			new[1].rotate(-1)
			new[1].glue(new[4])
			new[4].opposite()
			new[4].glue_as(b.copy().rotate(-1))
			self.Tetrahedra.remove(a.Tetrahedron)
			self.Tetrahedra.remove(b.Tetrahedron)
		else:
			return 0
		self.clear_and_rebuild()
		return 1









	








