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
		# Make sure all tets have at least the identity map in their Symmetries list, required by some functions.
		for tet in self.Tetrahedra:
			if len(tet.Symmetries) == 0:
				tet.Symmetries.append(Perm4((0,1,2,3)))
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
					sanity_check = 0
					while active:
						if len(completed_walks) > 0:
							shortest_completed_walk = completed_walks[0]
							for walk in completed_walks:
								if len(walk) < len(shortest_completed_walk):
									shortest_completed_walk = walk
						sanity_check = sanity_check + 1
						#print(sanity_check)
						if sanity_check > 1000:
							print('error building edge classes for starting arrow',first_arrow)
							break
						walk = active.pop()
						if len(completed_walks) == 0 or len(walk) < len(shortest_completed_walk):
							a = walk[-1]
							for sym in a.Tetrahedron.Symmetries:
								new_a = Arrow(sym.image(a.Edge),sym.image(a.Face),a.Tetrahedron)
								if new_a.next() != None:
									# now we want to check if new_a (which has been changed by .next()), or the image 
									# under some symmetry of new_a, is first_arrow. If so, then walk is a complete walk.
									complete = False
									if new_a.Tetrahedron == first_arrow.Tetrahedron:
										for other_sym in new_a.Tetrahedron.Symmetries:
											if other_sym.image(new_a.Edge) == first_arrow.Edge and other_sym.image(new_a.Face) == first_arrow.Face:
												complete = True
									if complete:
										completed_walks.append(walk)
									else:
										walk_branch = walk + [new_a]
										active.append(walk_branch)
					shortest_walk = completed_walks[0]
					for walk in completed_walks:
						if len(walk) < len(shortest_walk):
							shortest_walk = walk
					new_edge = Edge()
					z = ComplexSquareRootCombination.One()
					for a in shortest_walk:
						z = z*a.Tetrahedron.edge_params[a.Edge]
						new_edge.Corners.append(Corner(a.Tetrahedron,a.Edge))
						for perm in a.Tetrahedron.Symmetries:
							a.Tetrahedron.Class[perm.image(a.Edge)] = new_edge
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


	# set self.is_canonical to True or False
	def see_if_canonical(self):
		for tet1 in self.Tetrahedra:
			for face1 in TwoSubsimplices:
				if tet1.Neighbor[face1] != None:
					tet2 = tet1.Neighbor[face1]
					face2 = tet1.Gluing[face1].image(face1)
					if (tet1.tilt(comp(face1)) + tet2.tilt(comp(face2))).evaluate() > 0:
						self.is_canonical = False
						return
		# if it hasn't already returned, set it to True.
		self.is_canonical = True



	"""
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

	It could be that perm already does not descend to a well-defined map on the orbifold. We chech this with
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
	Simplification moves.
	"""

	"""
	2-3 move using arrows. It will automatically do a flat 2-3 move if necessary. There's an optional argument
	called build. By default, it is 1. But if it's set to 0, then the 2-3 move will be done without
	updating quotient edges, faces, or vertices, nor updating horotriangles.
	"""
	def arrow_two_to_three(self, two_subsimplex, tet, build = 1):
		#Third try. Have to be more careful with gluing to external faces when there are symmetries.
		if tet.Neighbor[two_subsimplex] is None:
			return 0
		flat_tets = []
		One = ComplexSquareRootCombination.One()
		self.PachnerPath.append([self.Tetrahedra,two_subsimplex,tet])
		if tet.face_glued_to_self(two_subsimplex):
			for one_subsimplex in OneSubsimplices:
				if tet.Gluing[two_subsimplex].image(one_subsimplex) == one_subsimplex:
					if is_subset(one_subsimplex,two_subsimplex):
						a = Arrow(one_subsimplex,two_subsimplex,tet)
						b = self.new_arrow()
						b.Tetrahedron.fill_edge_params(a.Tetrahedron.edge_params[a.Edge])
						break
			z = a.Tetrahedron.edge_params[a.simplex_south_head()]
			w = b.Tetrahedron.edge_params[b.simplex_north_tail()]
		else:
			a = Arrow(PickAnEdge[two_subsimplex], two_subsimplex, tet)
			b = a.glued()
		z = a.Tetrahedron.edge_params[a.simplex_south_head()]
		w = b.Tetrahedron.edge_params[b.simplex_north_tail()]
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
			self.see_if_canonical()
		return 1

	"""
	3-2 move. Returns 0 if a 3-2 move is not possible, otherwise does the move on self and returns 1.
	"""
	def three_to_two(self,edge):
		One = ComplexSquareRootCombination.One()
		a = Arrow(edge.Corners[0].Subsimplex,LeftFace[edge.Corners[0].Subsimplex],
				edge.Corners[0].Tetrahedron)
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
		elif edge.valence() == 3:
			for corner in edge.Corners:
				a = Arrow(corner.Subsimplex,LeftFace[corner.Subsimplex],
					corner.Tetrahedron)
				for sym in a.Tetrahedron.Symmetries:
					if sym.image(a.Edge) != a.Edge:
						return 0
					elif sym.tuple() != (0,1,2,3):
						face_glued_to_self = True
				if face_glued_to_self:
					break
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
			if a.glued().Tetrahedron is None:
				a.reverse()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				b.glue(a.glued())
				b.glue(a.glued())
			else:
				b.glue(a.glued())
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
			if a.Tetrahedron.face_glued_to_self(a.Face):
				b.glue(a.glued())
				b.glue(a.glued())
			else:
				b.glue(a.glued())
			a.reverse()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				c.glue(a.glued())
				c.glue(a.glued())
			else:
				c.glue(a.glued())
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
			if a.glued().Tetrahedron is None:
				a.reverse()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				b.glue(a.glued())
				b.glue(a.glued())
			else:
				b.glue(a.glued())
			a.opposite()
			if a.glued().Tetrahedron is None:
				a.reverse()
			a.next()
			b.rotate(-1)
			b.glue(a.opposite().glued())
			b.rotate(-1)
			b.glue(a.reverse().glued())
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
				if a.Tetrahedron.face_glued_to_self(a.Face):
					b.glue(a.glued())
					b.glue(a.glued())
				else:
					b.glue(a.glued())
				a.reverse()
				if a.Tetrahedron.face_glued_to_self(a.Face):
					c.glue(a.glued())
					c.glue(a.glued())
				else:
					c.glue(a.glued())
				b.rotate(-1)
				c.rotate(1)
				a.reverse().opposite().next()
		for corner in edge.Corners:
			if corner.Tetrahedron in self.Tetrahedra:
				self.Tetrahedra.remove(corner.Tetrahedron)
		for tet in self.Tetrahedra:
			tet.clear_Class()
			tet.horotriangles = {V0:None, V1:None, V2:None, V3:None}
		for i in range(len(self.Tetrahedra)):
			self.Tetrahedra[i].Index = i
		self.Edges = []
		self.Faces = []
		self.Vertices = []
		self.build_vertex_classes()
		self.build_edge_classes()
		self.add_cusp_cross_sections()
		for cusp in self.Vertices:
			self.normalize_cusp(cusp)
		self.see_if_canonical()
		return 1

	"""
	ATTENTION. I've now adjusted the 2-3 move above so that it also can do flat 2-3 moves. So the
	following function can maybe be deleted at some point.

	Flat 2-3 move. It could be that when you do a 2-3 move, one of the resulting tetrahedra is flat.
	Normally this is undesirable, so the 2-3 move should not be attempted. But if the flat tet corresponds
	to a pair of faces glued to themselves, you can just discard the flat tet. Return 0 if the flat 2-3 move
	doesn't make sense in this situation, otherwise do the move and return 1.
	"""
	def flat_two_to_three(self,two_subsimplex,tet):
		One = ComplexSquareRootCombination.One()
		for one_subsimplex in OneSubsimplices:
			if is_subset(one_subsimplex,two_subsimplex):
				z = tet.edge_params[one_subsimplex]
				w = tet.Neighbor[two_subsimplex].edge_params[tet.Gluing[two_subsimplex].image(one_subsimplex)]
				if (z*w).imag == SquareRootCombination.Zero():
					break
		if (z*w).imag != SquareRootCombination.Zero():
			return 0
		a = Arrow(one_subsimplex,two_subsimplex,tet)
		b = a.glued()
		a.reverse()
		if (not a.Tetrahedron.face_glued_to_self(a.Face) 
			or a.Tetrahedron.Gluing[a.Face].image(a.head()) != a.head()):
			return 0
		if (not b.Tetrahedron.face_glued_to_self(b.Face) 
			or b.Tetrahedron.Gluing[b.Face].image(b.head()) != b.head()):
			return 0
		#For now I'm not considering the case that two_subsimplex is glued to itself. Might want
		#to add that possibility later.
		z = b.Tetrahedron.edge_params[b.simplex_south_tail()]
		w = a.Tetrahedron.edge_params[a.simplex_axis()]
		c = self.new_arrow()
		d = self.new_arrow()
		c.Tetrahedron.fill_edge_params((One - z)*(One - w)/(z*w))
		d.Tetrahedron.fill_edge_params(z/(One - w))
		c.add_sym(c.copy())
		d.add_sym(d.copy())
		c.glue(d)
		d.glue(c)
		c.opposite()
		b.rotate(-1)
		if b.Tetrahedron.face_glued_to_self(b.Face):
			c.glue(b.glued())
			c.glue(b.glued())
		else:
			c.glue(b.glued())
		c.reverse()
		a.rotate(1)
		if a.Tetrahedron.face_glued_to_self(a.Face):
			c.glue(a.glued())
			c.glue(a.glued())
		else:
			c.glue(a.glued())
		d.opposite()
		b.rotate(-1)
		if b.Tetrahedron.face_glued_to_self(b.Face):
			d.glue(b.glued())
			d.glue(b.glued())
		else:
			d.glue(b.glued())
		d.reverse()
		a.rotate(1)
		if a.Tetrahedron.face_glued_to_self(a.Face):
			d.glue(a.glued())
			d.glue(a.glued())
		else:
			d.glue(a.glued())
		for tet in self.Tetrahedra:
			tet.clear_Class()
			tet.horotriangles = {V0:None, V1:None, V2:None, V3:None}
		self.Tetrahedra.remove(a.Tetrahedron)
		self.Tetrahedra.remove(b.Tetrahedron)
		for i in range(len(self.Tetrahedra)):
			self.Tetrahedra[i].Index = i
		self.Edges = []
		self.Faces = []
		self.Vertices = []
		self.build_vertex_classes()
		self.build_edge_classes()
		self.add_cusp_cross_sections()
		for cusp in self.Vertices:
			self.normalize_cusp(cusp)
		self.see_if_canonical()
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
	"""
	
	def three_to_six(self,two_subsimplex,tet):
		One = ComplexSquareRootCombination.One()
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
		if ((z*w*w).imag).evaluate() < 0 or (z*w*w).imag == SquareRootCombination.Zero():
			#In the latter case maybe a flat 3-6 move is possible, but let's not consider that currently.
			return 0
		flat_u0_u1_v1_w1 = False
		flat_u0_u1_v0_w1 = False
		z = tet.edge_params[a.simplex_south_head()]
		w = voisin.edge_params[b.simplex_south_tail()]
		if ((z*w).imag).evaluate() < 0:
			return 0
		if (z*w).imag == SquareRootCombination.Zero():
			flat_u0_u1_v1_w1 = True
		z = tet.edge_params[a.simplex_north_head()]
		w = voisin.edge_params[b.simplex_north_tail()]
		if ((z*w).imag).evaluate() < 0:
			return 0
		if (z*w).imag == SquareRootCombination.Zero():
			flat_u0_u1_v0_w1 = True
		c = self.new_arrow()
		c.glue(a)
		c.Tetrahedron.fill_edge_params(b.Tetrahedron.edge_params[b.Edge])
		if flat_u0_u1_v0_w1 and flat_u0_u1_v1_w1:
			return 0
		if flat_u0_u1_v0_w1:
			a.reverse()
			self.arrow_two_to_three(a.Face,a.Tetrahedron,0)
			b.reverse()
			new_b = b.copy()
			b.next().opposite().next()
			self.arrow_two_to_three(new_b.Face,new_b.Tetrahedron,0)
			new_b = b.copy()
			b.reverse().next()
			self.arrow_two_to_three(new_b.Face,new_b.Tetrahedron,0)
			c = b.copy().reverse().next()
			c.add_sym(c.copy().opposite())
			c.opposite().next()
			d = b.copy().reverse().rotate(-1).next()
			d.add_sym(d.copy().reverse())
			c.Tetrahedron.erase()
			self.Tetrahedra.remove(c.Tetrahedron)
		elif flat_u0_u1_v1_w1:
			self.arrow_two_to_three(a.Face,a.Tetrahedron,0)
			new_c = c.copy()
			c.next().opposite().next()
			self.arrow_two_to_three(new_c.Face,new_c.Tetrahedron,0)
			new_c = c.copy().reverse()
			c.next()
			self.arrow_two_to_three(new_c.Face,new_c.Tetrahedron,0)
			b = c.copy().next()
			b.add_sym(b.copy().opposite().reverse())
			b.opposite().next()
			d = c.copy().rotate(-1).next()
			d.add_sym(d.copy().opposite())
			c.Tetrahedron.erase()
			self.Tetrahedra.remove(c.Tetrahedron)
		else:
			x = tet.edge_params[a.simplex_north_tail()]
			y = voisin.edge_params[b.Edge]
			z = voisin.edge_params[b.simplex_south_tail()]
			w0 = One - One/(y - x*y + x)
			complex_abs_w0 = ComplexSquareRootCombination(abs(w0),SquareRootCombination.Zero())
			w1 = One - One/(x*z)
			complex_abs_w1 = ComplexSquareRootCombination(abs(w1),SquareRootCombination.Zero())
			if (w0/complex_abs_w0).real.evaluate() < (w1/complex_abs_w1).real.evaluate():
				#This is the case [u_0,w_1] "beneath" [u_1,w_0].
				self.arrow_two_to_three(two_subsimplex,tet,0)
				c.next()
				new_c = c.copy().reverse()
				c.opposite().next().reverse()
				self.arrow_two_to_three(new_c.Face,new_c.Tetrahedron,0)
				new_c = c.copy()
				c.reverse().next()
				self.arrow_two_to_three(new_c.Face,new_c.Tetrahedron,0)
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
				self.arrow_two_to_three(c.Face,c.Tetrahedron,0)
				b.reverse()
				new_b = b.copy()
				b.next().opposite().next()
				self.arrow_two_to_three(new_b.Face,new_b.Tetrahedron,0)
				new_b = b.copy().next()
				self.arrow_two_to_three(new_b.Face,new_b.Tetrahedron,0)
				#Now add symmetries.
				d = b.copy()
				b.next()
				b.add_sym(b.copy().opposite())
				e = b.copy().opposite().reverse().next()
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
				f.rotate(1)
				f.glue(b.glued())
				b.Tetrahedron.erase()
				self.Tetrahedra.remove(d.Tetrahedron)
				self.Tetrahedra.remove(b.Tetrahedron)
			else:
				#Otherwise the geodesics intersect.
				return 0
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
		self.see_if_canonical()
		return 1


	"""
	def three_to_six(self,two_subsimplex,tet):
		One = ComplexSquareRootCombination.One()
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
		z = tet.edge_params[one_subsimplex]
		w = voisin.edge_params[tet.Gluing[two_subsimplex].image(one_subsimplex)]
		if (z*w*w).imag == 0 or ((z*w*w).imag).evaluate() < 0:
			return 0
		for other_edge in OneSubsimplices:
			if is_subset(other_edge,two_subsimplex) and other_edge != one_subsimplex:
				z = tet.edge_params[other_edge]
				w = voisin.edge_params[tet.Gluing[two_subsimplex].image(other_edge)]
				if (z*w).imag == 0 or ((z*w).imag).evaluate() < 0:
					return 0
		a = Arrow(one_subsimplex,two_subsimplex,tet)
		b = a.glued()
		x = tet.edge_params[a.simplex_north_tail()]
		y = voisin.edge_params[b.Edge]
		z = voisin.edge_params[b.simplex_south_tail()]
		if (((One - One/(y - x*y + x))/abs(One - One/(y - x*y + x))).real.evaluate() < 
			((One - One/(x*z))/abs(One - One/(x*z))).real.evaluate()):
			new = self.new_arrows(4)
			#In terms of my write-up labeling, new[0] = [u0,u1,v1,w1], new[1] = [u0,u1,w0,w1],
			#new[2] = [u1,v0,w0,w1], new[3] = [v0,v1,w0,w1].
			new[0].Tetrahedron.fill_edge_params((x*(One - (One - x)*y + x))/((x - One)*y))
			new[1].Tetrahedron.fill_edge_params(((y - x*y + x)*(One - x*z))/(y - x*y + x - x*z))
			new[2].Tetrahedron.fill_edge_params((x*z - One)/(y - x*y + x - One))
			new[3].Tetrahedron.fill_edge_params((x - x*z)/(y - x*y + x - x*z))
			new[0].add_sym(new[0].copy())
			new[1].add_sym(new[1].copy())
			new[1].add_sym(new[1].copy().opposite())
			new[2].add_sym(new[2].copy())
			new[3].add_sym(new[3].copy())
			new[3].add_sym(new[3].copy().opposite())
			new[1].rotate(-1).reverse()
			new[0].glue(new[1])
			new[0].rotate(1)
			b.reverse().rotate(-1)
			if b.Tetrahedron.face_glued_to_self(b.Face):
				new[0].glue(b.glued())
				new[0].glue(b.glued())
			else:
				new[0].glue(b.glued())
			new[0].rotate(1)
			a.rotate(1)
			if a.glued().Tetrahedron is None:
				a.opposite().reverse()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				new[0].glue(a.glued())
				new[0].glue(a.glued())
			else:
				new[0].glue(a.glued())
			new[0].reverse()
			new[2].reverse()
			new[0].glue(new[2])
			b.opposite()
			if b.Tetrahedron.face_glued_to_self(b.Face):
				new[2].glue(b.glued())
				new[2].glue(b.glued())
			else:
				new[2].glue(b.glued())
			new[1].reverse().rotate(-1)
			new[2].reverse().opposite()
			new[1].glue(new[2])
			new[3].rotate(1)
			new[2].glue(new[3])
			new[3].reverse().rotate(-1)
			b.opposite().reverse()
			if b.Tetrahedron.face_glued_to_self(b.Face):
				new[3].glue(b.glued())
				new[3].glue(b.glued())
			else:
				new[3].glue(b.glued())
			#I think that's it for this case.


		elif (((One - One/(y - x*y + x))/abs(One - One/(y - x*y + x))).real.evaluate() > 
			((One - One/(x*z))/abs(One - One/(x*z))).real.evaluate()):

		else:
			return 0
	"""

	








