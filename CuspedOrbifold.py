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
		# Make sure all tets have at least the identity map in their Symmetries list, required by some functions.
		for tet in self.Tetrahedra:
			if len(tet.Symmetries) == 0:
				tet.Symmetries.append(Perm4((0,1,2,3)))
		self.Edges = []
		self.Faces = []
		self.Vertices = []
		self.DestSeq = None
		self.PachnerPath = []
		for T in self.Tetrahedra:
			T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
		self.build_vertex_classes()
		self.build_edge_classes()
		self.add_cusp_cross_sections()
		for cusp in self.Vertices:
			self.normalize_cusp(cusp)

	def info(self):
		tets = self.Tetrahedra
		print('number of tetrahedra is',len(tets))
		for i in range(len(tets)):
			print('gluing data for', tets[i], 'is')
			for j in range(4):
				print('face', j, 'of', tets[i], 'glued to', tets[i].Neighbor[TwoSubsimplices[j]])
				print('with gluing map',tets[i].Gluing[TwoSubsimplices[j]])
			print('symmetries of',tets[i],'are')
			for sym in tets[i].Symmetries:
				print(sym)
			print('edge_params of',tets[i],'are')
			print('edge 01',tets[i].edge_params[E01].real,'+',tets[i].edge_params[E01].imag,'* i')
			print('edge 02',tets[i].edge_params[E02].real,'+',tets[i].edge_params[E02].imag,'* i')
			print('edge 03',tets[i].edge_params[E03].real,'+',tets[i].edge_params[E03].imag,'* i')
			print('edge 12',tets[i].edge_params[E12].real,'+',tets[i].edge_params[E12].imag,'* i')
			print('edge 13',tets[i].edge_params[E13].real,'+',tets[i].edge_params[E13].imag,'* i')
			print('edge 23',tets[i].edge_params[E23].real,'+',tets[i].edge_params[E23].imag,'* i')

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
			for edge in OneSubsimplices:
				new_tet.edge_params[edge] = new_to_old[new_tet].edge_params[edge]
			for sym in new_to_old[new_tet].Symmetries:
				new_tet.Symmetries.append(sym)
			for face in TwoSubsimplices:
				if new_to_old[new_tet].Neighbor[face] != None:
					new_tet.attach(face,old_to_new[new_to_old[new_tet].Neighbor[face]],
						new_to_old[new_tet].Gluing[face].tuple())
		return CuspedOrbifold(new_tets)

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
			SquareRootCombination.One())
		active = [(tet0, vert0)]
		while active:
			tet0, vert0 = active.pop()
			for face0 in FacesAnticlockwiseAroundVertices[vert0]:
				# We need the following for orbifolds becuse some faces might be unglued.
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

		target_area = SquareRootCombination.SqrtThree()

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
	Here's my attempt at a new build_edge_classes function. Will delete the old one, which is below it,
	after enough testing. I'm re-doing it because it can be much simpler, and because I want to 
	allow for self to not be geometric (?). In that case we don't have edge_params and the edge locus
	order should already be given; it doesn't need to be figured out by this function.

	Two edges are in the same class if they're identified by a sequence of face gluings and
	symmetries. In that case, the two edges are assigned the same edge object (defined in edge.py).
	For example, if edge E01 of tet0 has the edge object EDGE, then tet0.Class[E01] == EDGE.
	One attribute of the edge class is "Corners". A Corner is an object with attributes Tetrahedron
	and Subsimplex. It might make sense for EDGE.Corners to be the list of all Corners such that
	Corner.Tetrahedron.Class[Corner.Subsimplex] == EDGE. But in fact we do things slightly
	differently. We walk around EDGE using arrows, and only the (one)Subsimplex-Tetrahedron pairs we
	encounter there get recorded into EDGE.Corners. Because of symmetries, this might not be all of them.
	Additionally, there could be repeats in this list, if an arrow and its reverse both apppear
	in the walk. We do it this way because we need to know the total dihedral angle traversed in a
	walk, theta. Then the integer n satisfying n*theta = 2*pi is the order of the local group fixing the
	edge, recorded as EDGE.LocusOrder == n.  
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
							print('Bad gluing data: could not construct edge link.')
							break
						newEdge.Corners.append(Corner(a.Tetrahedron, a.Edge))
						if a.true_next() is None:
							print('Hit boundary. Did not construct edge link.')
							break
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
					#We also take the product of all the edge_params in this walk.
					z = ComplexSquareRootCombination.One()
					for corner in newEdge.Corners:
						z = z*corner.Tetrahedron.edge_params[corner.Subsimplex]
						for sym in corner.Tetrahedron.Symmetries:
							corner.Tetrahedron.Class[sym.image(corner.Subsimplex)] = newEdge
					#Now we figure out the order of the local group fixing newEdge. This is
					#the smallest n such that z^n == 1.
					n = 0
					w = ComplexSquareRootCombination.One()
					while n < 1000:
						w = w*z
						n = n + 1
						# So w = z^n
						if w == ComplexSquareRootCombination.One():
							newEdge.LocusOrder = n
							break
					else:
						raise Exception('error getting edge locus order')
		for i in range(len(self.Edges)):
			self.Edges[i].Index = i

	# Checks if the triangulation is proto_canonical or not.
	def is_proto_canonical(self):
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
		copy = self.copy()
		# Make a list of flat tets.
		flat_tets = []
		for tet in copy.Tetrahedra:
			if tet.is_flat():
				flat_tets.append(tet)
		if len(flat_tets) == 0:
			return True
		for tet in flat_tets:
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
						raise Exception("Error when checking if proto_canonical")
					elif new_0.is_flat():
						if not (new_1.tilt(V3) + new_2.tilt(V2)).evaluate() > 0:
							# Then the edge is not concave, so the tet is not admissible.
							return False
					else:
						if not (new_1.tilt(V2) + new_0.tilt(V3)).evaluate() > 0:
							# Then the edge is not concave, so the tet is not admissible.
							return False
					# If we haven't returned False, then that flat tet is admissible. We should
					# reverse the special_four_to_four move before moving on to other flat tets. 
					copy.four_to_four(new_0.Class[E01])
					break
			else:
				# If we hit this, then this is a flat tet we weren't even able to do
				# a special_four_to_four move on, so it's definitely not admissible.
				return False
		return True


	"""
	I'M RE-DOING THIS IN THE NEW FILE simplicial_maps.py.

	Now I'm going to make the isometry_group and are_isometric functions. 

	The function valid_tet_to_tet(self,tet_a,tet_b,map_a_to_b) checks if mapping tet_a to tet_b with map_a_to_b
	respects symmetries and the locus orders of the edges.
    """

	def valid_tet_to_tet(self,tet_a,tet_b,map_a_to_b):
		for sym_a in tet_a.Symmetries:
			if (map_a_to_b*sym_a*inv(map_a_to_b)).tuple() not in [sym_b.tuple() for sym_b in tet_b.Symmetries]:
				return False
		for sym_b in tet_b.Symmetries:
			if (inv(map_a_to_b)*sym_b*map_a_to_b).tuple() not in [sym_a.tuple() for sym_a in tet_a.Symmetries]:
				return False
		for edge in OneSubsimplices:
			if tet_a.Class[edge].LocusOrder != tet_b.Class[map_a_to_b.image(edge)].LocusOrder:
				return False
		return True

	
	"""
	The next function tries to extend perm: tet0 --> teti to an isometry defined on every tet, given as a dictionary
	"isom" whose assignments are: if tetj is mapped to tetk via "phi" then isom[tetj] = (phi, tetk).

	It could be that perm already does not descend to a well-defined map on the orbifold. We check this with
	valid_tet_to_tet.

	If perm is valid, to extend perm we look at the neighbors of tet0 and extend the map to them in the unique way (unique
	up to post-composition with a symmetry), then check the maps on those tetrahedra are valid with valid_tet_to_tet,
	then add those tets to a queue. Then we pop the last element of the queue and see if its neighbors have 
	been mapped to something yet. If a neighbor hasn't, then we extend the map to it and add it to the queue. 
	If isom is already defined on a neighbor, then we leave it alone. After all neighbors are inspected, we pop 
	another tet off the queue and keep going until the queue is empty.

	At this point we have a map isom which is defined on every tet, but it might not descend to a well-defined
	map of the orbifold. In the last part of the function we check if it gives an actual isometry of the
	orbifold. If it does, we return isom. If not, we return None.
	"""
	def check_extends(self,perm,image_of_tet0):
		isom = {self.Tetrahedra[0]:(perm,image_of_tet0)}
		if self.valid_tet_to_tet(self.Tetrahedra[0],image_of_tet0,perm) is False:
			return
		for tet in self.Tetrahedra:
			if tet.Index > 0:
				isom[tet] = None
		active = [self.Tetrahedra[0]]
		while active:
			tet = active.pop()
			image_of_tet = isom[tet][1]
			for face in TwoSubsimplices:
				voisin = tet.Neighbor[face]
				if voisin != None:
					# if we've already defined what the isom should do to voisin, then we do nothing.
					# otherwise we now define what it does to voisin and add it to active.
					if isom[voisin] is None:
						phi = isom[tet][0]*inv(tet.Gluing[face])
						for sym in image_of_tet.Symmetries:
							if image_of_tet.Neighbor[sym.image(isom[tet][0].image(face))] != None:
								image_of_voisin = image_of_tet.Neighbor[sym.image(isom[tet][0].image(face))]
								isom[voisin] = (image_of_tet.Gluing[sym.image(isom[tet][0].image(face))]*sym*phi,image_of_voisin)
								if self.valid_tet_to_tet(voisin,image_of_voisin,isom[voisin][0]) is False:
									return
								break
						active.append(voisin)
		# check that isom respects the face gluings. By construction it respects some face gluings, but maybe not all.
		for tet1 in self.Tetrahedra:
			for face1 in TwoSubsimplices:
				if tet1.Neighbor[face1] != None:
					tet2,face2 = glued_to(tet1,face1)
					well_defined_on_face1_and_face2 = False
					for sym1 in isom[tet1][1].Symmetries:
						for sym2 in isom[tet2][1].Symmetries:
							image_of_face1 = (sym1*isom[tet1][0]).image(face1)
							image_of_face2 = (sym2*isom[tet2][0]).image(face2)
							if isom[tet1][1].Neighbor[image_of_face1] != None:
								if glued_to(isom[tet1][1],image_of_face1) == (isom[tet2][1],image_of_face2):
									if (isom[tet1][1].Gluing[image_of_face1]*sym1*isom[tet1][0]).tuple() == (sym2*isom[tet2][0]*tet1.Gluing[face1]).tuple():
										well_defined_on_face1_and_face2 = True
					if well_defined_on_face1_and_face2 is False:
						return
		# if we haven't returned None by now, then isom is a valid isometry.
		return isom


	"""
	Finally, the next function computes the isometry group of an orbifold using the above two functions.
	Really they're the combinatorial automorphisms of the triangulation. If the triangulation is canonical,
	they will be all isometries. It includes the orientation reversing isometries.

	The isometries are given as dictionaries as described in the above function's comments.
	"""

	def isometries(self):
		isometries = []
		seen_maps_of_tet0 = []
		for tet in self.Tetrahedra:
			for permutation in Perm4.S4():
				if (permutation.tuple(),tet) not in seen_maps_of_tet0:
					for sym in tet.Symmetries:
						seen_maps_of_tet0.append(((sym*permutation).tuple(),tet))
					isom = self.check_extends(permutation,tet)
					if isom != None:
						isometries.append(isom)
		return isometries


	"""
	If we just want the orientation preserving isometries.
	"""
	
	def isometriesOP(self):
		isometries = self.isometries()
		isometriesOP = [iso for iso in isometries if iso[self.Tetrahedra[0]][0].sign() == 0]
		return isometriesOP
	


	"""
	PACHNER MOVES.
	"""

	"""
	First, the function check_two_to_three(self,two_subsimplex,tet) checks if a 2-3 move through
	two_subsimplex of tet is possible. It is impossible if certain dihedral angles are too large
	or if tet has symmetries which don't preserve two_subsimplex.

	Next, the function two_to_three(self,two_subsimplex,tet) actually does the 2-3 move. I chose
	to implement it using arrows. It will automatically do a "flat" 2-3 move if possible. Also,
	there's an optional argument called "build". By default, it's 1. But if it's set to 0, then
	the 2-3 move will be done without updating quotient edges, faces, or vertices, nor updating
	horotriangles. This is useful when we do the 3-6 move further down.
	"""
	def check_two_to_three(self,two_subsimplex,tet):
		if tet.Neighbor[two_subsimplex] is None:
			return 0
		if tet.Neighbor[two_subsimplex] == tet and not tet.face_glued_to_self(two_subsimplex):
			return 0
		for one_subsimplex in OneSubsimplices:
			if is_subset(one_subsimplex,two_subsimplex):
				z = tet.edge_params[one_subsimplex]
				w = tet.Neighbor[two_subsimplex].edge_params[tet.Gluing[two_subsimplex].image(one_subsimplex)]
				if (z*w).imag.evaluate() < 0:
					return 0
				""""
				This part is used to make sure that, if a flat tet is created, we can cancel it
				right away. But maybe we won't be able to cancel it immediately, but will be able to
				cancel it later. So I don't want this part anymore, will delete it eventually.
				if (z*w).imag == SquareRootCombination.Zero():
					face = flip_face(one_subsimplex,two_subsimplex)
					if not tet.face_glued_to_self(face):
						return 0
					if tet.Gluing[face].image(one_subsimplex) != one_subsimplex:
						return 0
				"""
		for sym in tet.Symmetries:
			if sym.image(two_subsimplex) != two_subsimplex:
				return 0
		for sym in tet.Neighbor[two_subsimplex].Symmetries:
			if sym.image(tet.Gluing[two_subsimplex].image(two_subsimplex)) != tet.Gluing[two_subsimplex].image(two_subsimplex):
				return 0
		return 1

	

	"""
	Re-doing 2-3 move so it allows flat tets to be created.	
	"""
	def two_to_three(self, two_subsimplex, tet, build = 1):
		#Third try. Have to be more careful with gluing to external faces when there are symmetries.
		if tet.Neighbor[two_subsimplex] is None:
			return 0
		One = ComplexSquareRootCombination.One()
		if tet.face_glued_to_self(two_subsimplex):
			for one_subsimplex in OneSubsimplices:
				if tet.Gluing[two_subsimplex].image(one_subsimplex) == one_subsimplex:
					if is_subset(one_subsimplex,two_subsimplex):
						a = Arrow(one_subsimplex,two_subsimplex,tet)
						b = self.new_arrow()
						b.Tetrahedron.fill_edge_params(a.Tetrahedron.edge_params[a.Edge])
						break
		else:
			a = Arrow(PickAnEdge[two_subsimplex], two_subsimplex, tet)
			b = a.glued()
		z = a.Tetrahedron.edge_params[a.south_head()]
		w = b.Tetrahedron.edge_params[b.north_tail()]
		#Now we make the new tets. We consider each case separately. That means some redundant
		#lines are written, but I think it's easier to understand this way. 
		if tet.face_glued_to_self(two_subsimplex) and tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			new.Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			new.add_sym(new.copy().reverse())
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new.glue(new.copy().reverse())
			new.copy().reverse().glue(new)
			a.reverse()
			new.opposite()
			if a.glued().Tetrahedron is None:
				#Then must apply a rotation to see what the face of new is glued to.
				a.rotate(1)
				if a.glued().Tetrahedron is None:
					#Then when we rotate one more time, that face must be glued to something.
					a.rotate(1)
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new.glue(a.glued())
				new.glue(a.glued())
			else:
				new.glue(a.glued())
			#Now new is glued either to itself or to another tet which is not a.Tetrahedron
			#nor b.Tetrahedron. So all the data for new.Tetrahedron is assigned.
		elif tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			new.Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new.glue(new.copy())
			a.reverse()
			new.opposite()
			if a.glued().Tetrahedron is None:
				#Then must apply a rotation to see what the face of new is glued to.
				a.rotate(1)
				if a.glued().Tetrahedron is None:
					#Then when we rotate one more time, that face must be glued to something.
					a.rotate(1)
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new.glue(a.glued())
				new.glue(a.glued())
			else:
				new.glue(a.glued())
			new.reverse()
			if b.glued().Tetrahedron is None:
				#Then must apply a rotation to see what the face of new is glued to.
				b.rotate(1)
				if b.glued().Tetrahedron is None:
					#Then when we rotate one more time, that face must be glued to something.
					b.rotate(1)
			if b.Tetrahedron.face_glued_to_self(b.Face):
				new.glue(b.glued())
				new.glue(b.glued())
			else:
				new.glue(b.glued())
			#Now all the data for new is assigned.
		elif tet.face_glued_to_self(two_subsimplex):
			new = self.new_arrows(2)
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[0].Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			new[0].add_sym(new[0].copy().reverse())
			new[1].Tetrahedron.fill_edge_params(z/(One-w))
			new[0].glue(new[1])
			new[1].glue(new[1].copy().reverse())
			#Now we glue the external faces.
			a.reverse()
			new[0].opposite()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new[0].glue(a.glued())
				new[0].glue(a.glued())
			else:
				new[0].glue(a.glued())
			a.rotate(-1)
			new[1].opposite()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new[1].glue(a.glued())
				new[1].glue(a.glued())
			else:
				new[1].glue(a.glued())
			a.rotate(-1)
			new[1].reverse()
			#No faces of b.Tetrahedron are glued to anything, since it's a copy of a.Tetrahedron.
			#So we instead glue new[1] to what a is glued to.
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new[1].glue(a.glued())
				new[1].glue(a.glued())
			else:
				new[1].glue(a.glued())
		else:
			#Assuming there is no face rotation nor is the face glued to itself.
			#Then this is just a normal 2-3 move.
			new = self.new_arrows(3)
			new[0].Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			new[1].Tetrahedron.fill_edge_params(z/(One-w))
			new[2].Tetrahedron.fill_edge_params(w/(One-z))
			for i in range(3):
				new[i].glue(new[(i+1)%3])
			a.reverse()
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
				c.opposite()
				if a.Tetrahedron.face_glued_to_self(a.Face):
					c.glue(a.glued())
					c.glue(a.glued())
				else:
					c.glue(a.glued())
				c.reverse()
				if b.Tetrahedron.face_glued_to_self(b.Face):
					c.glue(b.glued())
					c.glue(b.glued())
				else:
					c.glue(b.glued())
				a.rotate(-1)
				b.rotate(1)		
		if tet in self.Tetrahedra:
			self.Tetrahedra.remove(tet)
		if b.Tetrahedron in self.Tetrahedra:
			self.Tetrahedra.remove(b.Tetrahedron)
		if build:
			self.clear_and_rebuild()
		return 1








	"""
	This one tries to remove flat tets right away.
	"""
	def older_two_to_three(self, two_subsimplex, tet, build = 1):
		#Third try. Have to be more careful with gluing to external faces when there are symmetries.
		if tet.Neighbor[two_subsimplex] is None:
			return 0
		flat_tets = []
		One = ComplexSquareRootCombination.One()
		if tet.face_glued_to_self(two_subsimplex):
			for one_subsimplex in OneSubsimplices:
				if tet.Gluing[two_subsimplex].image(one_subsimplex) == one_subsimplex:
					if is_subset(one_subsimplex,two_subsimplex):
						a = Arrow(one_subsimplex,two_subsimplex,tet)
						b = self.new_arrow()
						b.Tetrahedron.fill_edge_params(a.Tetrahedron.edge_params[a.Edge])
						break
		else:
			a = Arrow(PickAnEdge[two_subsimplex], two_subsimplex, tet)
			b = a.glued()
		z = a.Tetrahedron.edge_params[a.south_head()]
		w = b.Tetrahedron.edge_params[b.north_tail()]
		#Now we make the new tets. We consider each case separately. That means some redundant
		#lines are written, but I think it's easier to understand this way. 
		if tet.face_glued_to_self(two_subsimplex) and tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			new.Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			new.add_sym(new.copy().reverse())
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new.glue(new.copy().reverse())
			new.copy().reverse().glue(new)
			a.reverse()
			new.opposite()
			if a.glued().Tetrahedron is None:
				#Then must apply a rotation to see what the face of new is glued to.
				a.rotate(1)
				if a.glued().Tetrahedron is None:
					#Then when we rotate one more time, that face must be glued to something.
					a.rotate(1)
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new.glue(a.glued())
				new.glue(a.glued())
			else:
				new.glue(a.glued())
			#Now new is glued either to itself or to another tet which is not a.Tetrahedron
			#nor b.Tetrahedron. So all the data for new.Tetrahedron is assigned.
		elif tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			new.Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new.glue(new.copy())
			a.reverse()
			new.opposite()
			if a.glued().Tetrahedron is None:
				#Then must apply a rotation to see what the face of new is glued to.
				a.rotate(1)
				if a.glued().Tetrahedron is None:
					#Then when we rotate one more time, that face must be glued to something.
					a.rotate(1)
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new.glue(a.glued())
				new.glue(a.glued())
			else:
				new.glue(a.glued())
			new.reverse()
			if b.glued().Tetrahedron is None:
				#Then must apply a rotation to see what the face of new is glued to.
				b.rotate(1)
				if b.glued().Tetrahedron is None:
					#Then when we rotate one more time, that face must be glued to something.
					b.rotate(1)
			if b.Tetrahedron.face_glued_to_self(b.Face):
				new.glue(b.glued())
				new.glue(b.glued())
			else:
				new.glue(b.glued())
			#Now all the data for new is assigned.
		elif tet.face_glued_to_self(two_subsimplex):
			new = self.new_arrows(2)
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[0].Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			new[0].add_sym(new[0].copy().reverse())
			new[1].Tetrahedron.fill_edge_params(z/(One-w))
			new[0].glue(new[1])
			new[1].glue(new[1].copy().reverse())
			for i in range(2):
				if new[i].Tetrahedron.edge_params[new[i].Edge].imag == SquareRootCombination.Zero():
					flat_tets.append(new[i].copy())
			#Now we glue the external faces.
			a.reverse()
			new[0].opposite()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new[0].glue(a.glued())
				new[0].glue(a.glued())
			else:
				new[0].glue(a.glued())
			a.rotate(-1)
			new[1].opposite()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new[1].glue(a.glued())
				new[1].glue(a.glued())
			else:
				new[1].glue(a.glued())
			a.rotate(-1)
			new[1].reverse()
			#No faces of b.Tetrahedron are glued to anything, since it's a copy of a.Tetrahedron.
			#So we instead glue new[1] to what a is glued to.
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new[1].glue(a.glued())
				new[1].glue(a.glued())
			else:
				new[1].glue(a.glued())
		else:
			#Assuming there is no face rotation nor is the face glued to itself.
			#Then this is just a normal 2-3 move.
			new = self.new_arrows(3)
			new[0].Tetrahedron.fill_edge_params((One-z)*(One-w)/(z*w))
			new[1].Tetrahedron.fill_edge_params(z/(One-w))
			new[2].Tetrahedron.fill_edge_params(w/(One-z))
			for i in range(3):
				if new[i].Tetrahedron.edge_params[new[i].Edge].imag == SquareRootCombination.Zero():
					flat_tets.append(new[i].copy())
			for i in range(3):
				new[i].glue(new[(i+1)%3])
			a.reverse()
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
				c.opposite()
				if a.Tetrahedron.face_glued_to_self(a.Face):
					c.glue(a.glued())
					c.glue(a.glued())
				else:
					c.glue(a.glued())
				c.reverse()
				if b.Tetrahedron.face_glued_to_self(b.Face):
					c.glue(b.glued())
					c.glue(b.glued())
				else:
					c.glue(b.glued())
				a.rotate(-1)
				b.rotate(1)
		#Now, if there's a flat tet, we try to remove it. You can do this if its two external faces
		#are glued to themselves in the correct way. If at least the one in tet is glued to itself
		#properly, the following will go through.
		if flat_tets:
			print('doing 2-3 with a flat tet')
			a = flat_tets[0].copy()
			a.opposite()
			if (a.Tetrahedron.face_glued_to_self(a.Face) and 
				a.Tetrahedron.Gluing[a.Face].image(a.Edge) == a.Edge):
				if tet.face_glued_to_self(two_subsimplex):
					#Then it must be that a.Tetrahedron is the one which got the nontrivial symmetry.
					#In that case, when we remove it, there's one new tet left. Its face which was attached
					#to a.Tetrahedron must now be glued to itself.
					a.opposite()
					if a.glued().Tetrahedron is None:
						a.reverse().next().reverse()
					else:
						a.next().reverse()
					a.glue(a.copy().reverse())
				else:
					a.opposite().next().reverse()
					a.glue(a.copy().next().next())
				self.Tetrahedra.remove(flat_tets[0].Tetrahedron)
			else:
				#Then the face isn't glued to itself properly, so a flat 2-3 move isn't possible. If
				#check_two_to_three was passed, this shouldn't happen. Note that self has been altered
				#in the attempt.
				return 0
		if tet in self.Tetrahedra:
			self.Tetrahedra.remove(tet)
		if b.Tetrahedron in self.Tetrahedra:
			self.Tetrahedra.remove(b.Tetrahedron)
		if build:
			self.Edges = []
			self.Faces = []
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
			#Because b.Tetrahedron has rotations, we are done specifying its face pairings.
			#Except if the face we just glued is glued to itself. Then the other faces will
			#also be glued to themselves, and my convention is that I should do those gluings,
			#even though they're implied by the symmetries. 
			if b.Tetrahedron.face_glued_to_self(b.Face):
				perm = b.Tetrahedron.Gluing[b.Face]
				for sym in b.Tetrahedron.Symmetries:
					if sym.tuple() != (0,1,2,3):
						b.Tetrahedron.attach(sym.image(b.Face),b.Tetrahedron,(sym*perm*inv(sym)).tuple())
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
			#Now if either of these faces was glued to itself, we do that for the others determined
			#by the symmetries.
			if b.Tetrahedron.face_glued_to_self(b.Face):
				perm = b.Tetrahedron.Gluing[b.Face]
				for sym in b.Tetrahedron.Symmetries:
					if sym.tuple() != (0,1,2,3):
						b.Tetrahedron.attach(sym.image(b.Face),b.Tetrahedron,(sym*perm*inv(sym)).tuple())
			if c.Tetrahedron.face_glued_to_self(c.Face):
				perm = c.Tetrahedron.Gluing[c.Face]
				for sym in c.Tetrahedron.Symmetries:
					if sym.tuple() != (0,1,2,3):
						c.Tetrahedron.attach(sym.image(c.Face),c.Tetrahedron,(sym*perm*inv(sym)).tuple())
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
			a.opposite()
			if a.glued().Tetrahedron is None:
				a.reverse()
			a.next()
			b.rotate(-1)
			b.glue_as(a.opposite())
			b.rotate(-1)
			b.glue_as(a.reverse())
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
			for i in range(3):
				a.opposite()
				b.glue_as(a)
				a.reverse()
				c.glue_as(a)
				b.rotate(-1)
				c.rotate(1)
				a.reverse().opposite().next()
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
		for sym in tet.Symmetries:
			if sym.tuple() != (0,1,2,3):
				break
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
		"""
		if ((z*w*w).imag).evaluate() < 0 or (z*w*w).imag == SquareRootCombination.Zero():
			#In the latter case maybe a flat 3-6 move is possible, but let's not consider that currently.
			return 0
		"""
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
			x = tet.edge_params[a.north_tail()]
			y = voisin.edge_params[b.Edge]
			z = voisin.edge_params[b.south_tail()]
			w0 = One - One/(y - x*y + x)
			complex_abs_w0 = ComplexSquareRootCombination(abs(w0),SquareRootCombination.Zero())
			w1 = One - One/(x*z)
			complex_abs_w1 = ComplexSquareRootCombination(abs(w1),SquareRootCombination.Zero())
			"""
			if (w0/complex_abs_w0).real == (w1/complex_abs_w1).real:
				#In this case the geodesics [u0,w1] and [u1,w0] intersect. You might think that 
				#this could be allowed.
				#It would create a flat tet, but we create flat tets sometimes, hoping to cancel
				#them later. But I think it really doesn't make sense here for a different reason.
				#When you create the flat tet, its nontrivial symmetry is incompatible with the
				#way it's flat, i.e. where the edges with angle pi are. So we don't allow this case.
				print("[u0,w1] and [u1,w0] intersect")
				#return 0
			"""
			c = self.new_arrow()
			c.glue(a)
			c.Tetrahedron.fill_edge_params(b.Tetrahedron.edge_params[b.Edge])
			if ((w0/complex_abs_w0).real.evaluate() < (w1/complex_abs_w1).real.evaluate() or
				(w0/complex_abs_w0).real == (w1/complex_abs_w1).real):
				#This is the case [u_0,w_1] "beneath" [u_1,w_0].
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
					if e.Tetrahedron.face_glued_to_self(e.Face):
						d.glue(e.glued())
						d.glue(e.glued())
					else:
						d.glue(e.glued())
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
					if d.Tetrahedron.face_glued_to_self(d.Face):
						e.glue(d.glued())
						e.glue(d.glued())
					else:
						e.glue(d.glued())
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
			#Now the arrows are as in the picture.
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
			#Now the arrows are as in the picture.
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









	








