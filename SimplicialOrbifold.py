from tetrahedron import*
from arrow import*
from edge import*
from vertex import*
from face import*
from simplex import*
from corner import*
from perm4 import*
from CanonizeInfo import*


class SimplicialOrbifold:
	def __init__(self,tets_list):
		self.Tetrahedra = tets_list
		for i in range(len(self.Tetrahedra)):
			self.Tetrahedra[i].Index = i
		self.Edges = []
		self.Faces = []
		self.Vertices = []
		self.build_vertex_classes()
		self.build_edge_classes()
		self.build_face_classes()

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
			print('edge labels of',tets[i],'are')
			print('edge 01: ',tets[i].edge_labels[E01])
			print('edge 02: ',tets[i].edge_labels[E02])
			print('edge 03: ',tets[i].edge_labels[E03])
			print('edge 12: ',tets[i].edge_labels[E12])
			print('edge 13: ',tets[i].edge_labels[E13])
			print('edge 23: ',tets[i].edge_labels[E23])

	def add_tet(self, tet):
		self.Tetrahedra.append(tet)

	# Add one new tetrahedron and return one of its arrows.
	def new_arrow(self):
		tet = Tetrahedron()
		self.add_tet(tet)
		return Arrow(E01,F3,tet)

	# Or, add a whole bunch of them.
	def new_arrows(self,n):
		return [self.new_arrow() for i in range(n)]

	def clear_and_rebuild(self):
		self.Edges = []
		self.Vertices = []
		for tet in self.Tetrahedra:
			tet.clear_Class()
		for i in range(len(self.Tetrahedra)):
			self.Tetrahedra[i].Index = i
		self.build_vertex_classes()
		self.build_edge_classes()

	def build_face_classes(self):
		for tet in self.Tetrahedra:
			for two_subsimplex in TwoSubsimplices:
				if tet.Class[two_subsimplex] == None:
					nbr, perm = tet.true_glued(two_subsimplex)
					nbr_two_subsimplex = perm.image(two_subsimplex)
					newFace = Face()
					self.Faces.append(newFace)
					for sym in tet.Symmetries:
						image = sym.image(two_subsimplex)
						if tet.Class[image] == None:
							tet.Class[image] = newFace
							newFace.Corners.append(Corner(tet,image))
					for sym in nbr.Symmetries:
						image = sym.image(nbr_two_subsimplex)
						if nbr.Class[image] == None:
							nbr.Class[image] = newFace
							newFace.Corners.append(Corner(nbr,image))
		for i in range(len(self.Faces)):
			self.Faces[i].Index = i

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
	We create edge classes just like we do for a CuspedOrbifold. The only difference is
	that we don't have to figure out the LocusOrder at the end using the geometry.

	Note that edge.Corners will not give all 2-simplices belonging to the edge class. It only
	contains the 2-simplices encountered on a complete walk around the edge. Some other 2-simplices
	which are in the class due to symmetries will not be in it. This is not how we do things for
	the vertex or face classes; there, all simplices in the class are in Corners. The reason is that
	edge.valence(), which is the length of edge.Corners, needs to equal the "length" of the walk
	around the edge, i.e. the number of distinct ARROWS around the edge (distint up to symmetries). 

	Also, we should allow for boundary faces probably. It currently doesn't do that, but I might
	change it. We should allow boundary faces because we might want to get the isomorphism signature
	of a simplicial orbifold with boundary, like if we're trying to build all triangulations
	up to a certain amount of tetrahedra and we want to avoid paths we've already gone down.
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
						if sanity_check > 6*len(self.Tetrahedra):
							print('Bad gluing data: could not construct edge link.')
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
					for corner in newEdge.Corners:
						for sym in corner.Tetrahedron.Symmetries:
							corner.Tetrahedron.Class[sym.image(corner.Subsimplex)] = newEdge
					#Now let's assign the LocusOrder of newEdge.
					newEdge.LocusOrder = corner.Tetrahedron.edge_labels[corner.Subsimplex]
		for i in range(len(self.Edges)):
			self.Edges[i].Index = i

	def change_edge_labels(self,edge,new_label):
		for corner in edge.Corners:
			tet = corner.Tetrahedron
			one_subsimplex = corner.Subsimplex
			for sym in tet.Symmetries:
				tet.edge_labels[sym.image(one_subsimplex)] = new_label
		edge.LocusOrder = new_label


	"""
	PACHNER MOVES.

	Here we define functions which do the 2-3, 3-2, 1-4, 4-1, 0-2, and 2-0 orbifold Pachner moves.
	They're similar to how they're written in the CuspedOrbifold class, except there's
	no geometry. And of course there are no 1-4 or 4-1 moves in the Cusped Orbifold class since
	they add or remove a non-ideal vertex.
	"""

	"""
	2-3 move.
	"""
	def check_two_to_three(self,two_subsimplex,tet):
		# Check if tet or its neighbor have symmetries making it impossible to do a 2-3 move.
		if len(tet.Symmetries) != 1 and len(tet.Symmetries) != 3:
			return False
		nbr = tet.Neighbor[two_subsimplex]
		if len(tet.Symmetries) != len(nbr.Symmetries):
			return False
		perm = tet.Gluing[two_subsimplex]
		nbr_face = perm.image(two_subsimplex)
		if len(tet.Symmetries) == 3:
			for sym in tet.Symmetries:
				if sym.image(two_subsimplex) != two_subsimplex:
					return False
			for sym in nbr.Symmetries:
				if sym.image(nbr_face) != nbr_face:
					return False
		return True

	def two_to_three(self, two_subsimplex, tet, build = 1):
		if tet.Neighbor[two_subsimplex] is None:
			return 0
		if tet.face_glued_to_self(two_subsimplex):
			for one_subsimplex in OneSubsimplices:
				if tet.Gluing[two_subsimplex].image(one_subsimplex) == one_subsimplex:
					if is_subset(one_subsimplex,two_subsimplex):
						a = Arrow(one_subsimplex,two_subsimplex,tet)
						break
		else:
			a = Arrow(PickAnEdge[two_subsimplex], two_subsimplex, tet)
			b = a.glued()
		# Now we make the new tets. We consider each case separately. It could be done more
		# quickly, but I think it's easier to understand this way. 
		if tet.face_glued_to_self(two_subsimplex) and tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new.add_sym(new.copy().reverse())
			new.glue(new.copy().reverse())
			a.reverse()
			new.opposite()
			new.true_glue_as(a)
			# That takes care of face gluings and symmetries. Now we need to
			# give the new tet the correct edge labels.
			new.Tetrahedron.edge_labels = {
				new.equator(): 3, 
				new.axis(): tet.edge_labels[a.axis()],
				new.north_head(): tet.edge_labels[a.north_head()], 
				new.north_tail(): tet.edge_labels[a.south_head()],
				new.south_head(): tet.edge_labels[a.south_head()],
				new.south_tail(): tet.edge_labels[a.north_head()]
				}
		elif tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new.glue(new.copy())
			a.reverse()
			new.opposite().true_glue_as(a)
			new.reverse().true_glue_as(b)
			# Now we set the edge labels.
			new.Tetrahedron.edge_labels = {
				new.equator(): 3,
				new.axis(): tet.edge_labels[a.axis()],
				new.north_head(): b.Tetrahedron.edge_labels[b.north_head()],
				new.north_tail(): tet.edge_labels[a.south_head()],
				new.south_head(): b.Tetrahedron.edge_labels[b.south_head()],
				new.south_tail(): tet.edge_labels[a.north_head()] 
				}
			# Update canonize info.
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
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[0].add_sym(new[0].copy().reverse())
			new[0].glue(new[1])
			new[1].glue(new[1].copy().reverse())
			# Now we glue the external faces.
			a.reverse()
			new[0].opposite().glue_as(a)
			a.rotate(-1)
			new[1].opposite().glue_as(a)
			a.rotate(-1)
			new[1].reverse().glue_as(a)
			# Now the edge labels.
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
			# Assuming there is no face rotation nor is the face glued to itself.
			# Then this is just a normal 2-3 move.
			new = self.new_arrows(3)
			for i in range(3):
				new[i].glue(new[(i+1)%3])
			a.reverse()
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
				c.opposite()
				c.glue_as(a)
				c.reverse()
				c.glue_as(b)
				a.rotate(-1)
				b.rotate(1)
			# Now the edge labels.
			for c in new:
				c.Tetrahedron.edge_labels = {
					c.equator(): 1,
					c.axis(): tet.edge_labels[a.axis()],
					c.north_head(): b.Tetrahedron.edge_labels[b.north_head()],
					c.north_tail(): tet.edge_labels[a.south_head()],
					c.south_head(): b.Tetrahedron.edge_labels[b.south_head()],
					c.south_tail(): tet.edge_labels[a.north_head()]
					}
				a.rotate(-1)
				b.rotate(1)
			# Update the canonize_info.
			if tet.canonize_info is not None:	
				for c in new:
					c.Tetrahedron.canonize_info = CanonizeInfo()
					c.Tetrahedron.canonize_info.part_of_coned_cell = True
					c.Tetrahedron.canonize_info.is_flat = False
					c.Tetrahedron.canonize_info.face_status = {
						c.north_face(): 2,
						c.south_face(): 2,
						c.east_face(): b.Tetrahedron.canonize_info.face_status[b.east_face()],
						c.west_face(): a.Tetrahedron.canonize_info.face_status[a.east_face()]
						}
					a.rotate(-1)
					b.rotate(1)
		self.Tetrahedra.remove(tet)
		if not tet.face_glued_to_self(two_subsimplex):
			# Then we need to remove b.Tetrahedron as well. If the face was glued to itself,
			# then there would be no b.
			self.Tetrahedra.remove(b.Tetrahedron)
		if build:
			self.clear_and_rebuild()
		return 1



	"""
	3-2 move.
	"""
	def three_to_two(self,edge,build = 1):
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
			a.opposite()
			b.true_glue_as(a)
			#Now the edge labels.
			b.Tetrahedron.edge_labels = {
					b.equator(): a.Tetrahedron.edge_labels[a.north_head()],
					b.axis(): a.Tetrahedron.edge_labels[a.axis()],
					b.north_head(): a.Tetrahedron.edge_labels[a.north_head()],
					b.north_tail(): a.Tetrahedron.edge_labels[a.axis()],
					b.south_head(): a.Tetrahedron.edge_labels[a.north_head()],
					b.south_tail(): a.Tetrahedron.edge_labels[a.axis()]
					}
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
			a.opposite()
			b.true_glue_as(a)
			a.reverse()
			c.true_glue_as(a)
			a.reverse()
			for d in [b,c]:
				d.Tetrahedron.edge_labels = {
					d.equator(): a.Tetrahedron.edge_labels[a.north_head()],
					d.axis(): a.Tetrahedron.edge_labels[a.axis()],
					d.north_head(): a.Tetrahedron.edge_labels[a.north_head()],
					d.north_tail(): a.Tetrahedron.edge_labels[a.axis()],
					d.south_head(): a.Tetrahedron.edge_labels[a.north_head()],
					d.south_tail(): a.Tetrahedron.edge_labels[a.axis()]
					}
				a.reverse() 
		elif face_glued_to_self:
			b = self.new_arrow()
			b.add_sym(b.copy())
			b.glue(b.copy().reverse())
			b.reverse()
			a.opposite()
			b.true_glue_as(a)
			c = a.copy().opposite()
			if c.glued().Tetrahedron is None:
				c.reverse()
			c.next()
			b.rotate(-1)
			b.glue_as(c.opposite())
			b.rotate(-1)
			b.glue_as(c.reverse())
			b.Tetrahedron.edge_labels = {
				b.equator(): c.Tetrahedron.edge_labels[c.south_tail()],
				b.axis(): a.Tetrahedron.edge_labels[a.axis()],
				b.north_head(): a.Tetrahedron.edge_labels[a.north_head()],
				b.north_tail(): c.Tetrahedron.edge_labels[c.axis()],
				b.south_head(): a.Tetrahedron.edge_labels[a.south_head()],
				b.south_tail(): c.Tetrahedron.edge_labels[c.axis()]
				}
			# This is the only case of the 3-2 move which could be used in canonize part 2.
			# For that situation, update canonize info.
			if a.Tetrahedron.canonize_info is not None:
				b.Tetrahedron.canonize_info = CanonizeInfo()
				b.Tetrahedron.canonize_info.part_of_coned_cell = True
				b.Tetrahedron.canonize_info.is_flat = False
				b.Tetrahedron.canonize_info.face_status = {
					b.north_face(): c.Tetrahedron.canonize_info.face_status[c.west_face()],
					b.south_face(): c.Tetrahedron.canonize_info.face_status[c.east_face()],
					b.east_face(): a.Tetrahedron.canonize_info.face_status[a.east_face()],
					b.west_face(): 2
					}
		else:
			b = self.new_arrow()
			c = self.new_arrow()
			b.add_sym(b.copy())
			c.add_sym(c.copy())
			b.glue(c)
			b.reverse()
			for i in range(3):
				a.opposite()
				b.glue_as(a)
				a.reverse()
				c.glue_as(a)
				a.reverse()
				b.Tetrahedron.edge_labels[b.axis()] = a.Tetrahedron.edge_labels[a.axis()]
				b.Tetrahedron.edge_labels[b.north_head()] = a.Tetrahedron.edge_labels[a.north_head()]
				c.Tetrahedron.edge_labels[c.axis()] = a.Tetrahedron.edge_labels[a.axis()]
				c.Tetrahedron.edge_labels[c.south_head()] = a.Tetrahedron.edge_labels[a.north_tail()]
				b.rotate(-1)
				c.rotate(1)
				a.opposite().next()
		for corner in edge.Corners:
			if corner.Tetrahedron in self.Tetrahedra:
				self.Tetrahedra.remove(corner.Tetrahedron)
		if build:
			self.clear_and_rebuild()
		return 1




	"""
	1-4 move. There's a different case for each kind of symmetry group (up to iso).

	The tetrahedron to be subdvided is tet. We modify self accordingly then return None.
	"""
	def one_to_four(self,tet):
		if len(tet.Symmetries) == 1:
			# Then there is only the trivial symmetry.
			# See the diagram in the notes to understand where the arrows are.
			new = self.new_arrows(4)
			for i in range(4):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			# First we make all the gluings.
			old = Arrow(E02,F3,tet)
			for i in range(3):
				new[i].glue(new[(i+1)%3])
				new[3].glue(new[i].copy().opposite())
				new[3].rotate(-1)
				new[i].copy().opposite().glue_as(old)
				old.rotate(-1)
			new[3].copy().reverse().glue_as(old.copy().reverse())
			# Now we set the edge labels.
			for i in range(3):
				new[i].Tetrahedron.edge_labels = {
					new[i].equator(): tet.edge_labels[old.axis()],
					new[i].axis(): 1,
					new[i].north_head(): tet.edge_labels[old.north_head()],
					new[i].north_tail(): tet.edge_labels[old.south_head()],
					new[i].south_head(): 1,
					new[i].south_tail(): 1 
					} 
				old.rotate(-1)
			new[3].Tetrahedron.edge_labels = {
				new[3].equator(): 1,
				new[3].axis(): tet.edge_labels[old.axis()],
				new[3].north_head(): 1,
				new[3].north_tail(): tet.edge_labels[old.north_tail()],
				new[3].south_head(): 1,
				new[3].south_tail(): tet.edge_labels[old.south_tail()] 
				}
			# Now update the canonize info.
			if tet.canonize_info is not None:
				for i in range(4):
					new[i].Tetrahedron.canonize_info = CanonizeInfo()
					new[i].Tetrahedron.canonize_info.part_of_coned_cell = True
					new[i].Tetrahedron.canonize_info.is_flat = False
				for i in range(3):
					new[i].Tetrahedron.canonize_info.face_status = {
						new[i].north_face(): tet.canonize_info.face_status[old.east_face()],
						new[i].south_face(): 2,
						new[i].east_face(): 2,
						new[i].west_face(): 2
						}
					old.rotate(-1)
				new[3].Tetrahedron.canonize_info.face_status = {
					new[3].north_face(): 2,
					new[3].south_face(): 2,
					new[3].east_face(): 2,
					new[3].west_face(): tet.canonize_info.face_status[old.west_face()]
				}
		if len(tet.Symmetries) == 2:
			new = self.new_arrows(2)
			for i in range(2):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			sym = tet.nontrivial_sym()
			for one_subsimplex in OneSubsimplices:
				if sym.image(one_subsimplex) == one_subsimplex:
					old = Arrow(one_subsimplex,RightFace[one_subsimplex],tet)
					break
			# Face gluings.
			new[0].glue(new[1])
			new[0].copy().opposite().reverse().glue(new[0].copy().opposite())
			new[0].copy().reverse().glue(new[1].copy().opposite().rotate(-1))
			new[0].copy().opposite().true_glue_as(old)
			new[1].glue(new[1].copy().reverse().rotate(-1))
			new[1].copy().opposite().true_glue_as(old.copy().rotate(-1))
			# Edge labels.
			for i in range(2):
				new[i].Tetrahedron.edge_labels = {
					new[i].equator(): tet.edge_labels[old.axis()],
					new[i].axis(): 1,
					new[i].north_head(): tet.edge_labels[old.north_head()],
					new[i].north_tail(): tet.edge_labels[old.south_head()],
					new[i].south_head(): 1,
					new[i].south_tail(): 1 
					} 
				old.rotate(-1)
			# Restore old to its position in the diagram.
			old.rotate(-1)
			# Update the canonize info.
			if tet.canonize_info is not None:
				for i in range(2):
					new[i].Tetrahedron.canonize_info = CanonizeInfo()
					new[i].Tetrahedron.canonize_info.part_of_coned_cell = True
					new[i].Tetrahedron.canonize_info.is_flat = False
				for i in range(2):
					new[i].Tetrahedron.canonize_info.face_status = {
						new[i].north_face(): tet.true_face_status(old.east_face()),
						new[i].south_face(): 2,
						new[i].east_face(): 2,
						new[i].west_face(): 2
						}
					old.rotate(-1)
		if len(tet.Symmetries) == 3:
			new = self.new_arrows(2)
			for i in range(2):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[0].add_sym(new[0].copy().rotate(1).reverse())
			new[0].add_sym(new[0].copy().reverse().rotate(-1))
			sym = tet.nontrivial_sym()
			for face in TwoSubsimplices:
				if sym.image(face) == face:
					break
			for edge in OneSubsimplices:
				if is_subset(edge,face):
					old = Arrow(edge,face,tet)
					break
			# Face gluings.
			new[0].glue(new[1])
			new[0].copy().opposite().glue_as(old)
			new[1].copy().reverse().rotate(-1).glue(new[1].copy().reverse().rotate(-1))
			new[1].copy().opposite().true_glue_as(old.copy().rotate(-1))
			# Edge labels.
			for i in range(2):
				new[i].Tetrahedron.edge_labels = {
					new[i].equator(): tet.edge_labels[old.axis()],
					new[i].axis(): 1,
					new[i].north_head(): tet.edge_labels[old.north_head()],
					new[i].north_tail(): tet.edge_labels[old.south_head()],
					new[i].south_head(): 1,
					new[i].south_tail(): 1 
					} 
				old.rotate(-1)
			# Actually, one of the edges of new[1].Tetrahedron gets labeled 3 because of the symmetry.
			new[1].Tetrahedron.edge_labels[new[1].south_head()] = 3
			# Return old to its diagram position.
			old.rotate(-1)
			# Canonize info.
			if tet.canonize_info is not None:
				for i in range(2):
					new[i].Tetrahedron.canonize_info = CanonizeInfo()
					new[i].Tetrahedron.canonize_info.part_of_coned_cell = True
					new[i].Tetrahedron.canonize_info.is_flat = False
				for i in range(2):
					new[i].Tetrahedron.canonize_info.face_status = {
						new[i].north_face(): tet.true_face_status(old.east_face()),
						new[i].south_face(): 2,
						new[i].east_face(): 2,
						new[i].west_face(): 2
						}
					old.rotate(-1)
		if len(tet.Symmetries) == 4:
			new = self.new_arrow()
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			old = Arrow(E02,F3,tet)
			# Face gluings.
			for i in range(3):
				new.glue(new.copy().reverse())
				new.rotate(1)
			new.reverse().true_glue_as(old)
			# Edge labels.
			# Move new to a different position so we can assign edge labels in the
			# same way as the other cases.
			new.opposite()
			new.Tetrahedron.edge_labels = {
				new.equator(): tet.edge_labels[old.axis()],
				new.axis(): 1,
				new.north_head(): tet.edge_labels[old.north_head()],
				new.north_tail(): tet.edge_labels[old.south_head()],
				new.south_head(): 1,
				new.south_tail(): 1 
				}
			# Canonize info.
			if tet.canonize_info is not None:
				new.Tetrahedron.canonize_info = CanonizeInfo()
				new.Tetrahedron.canonize_info.part_of_coned_cell = True
				new.Tetrahedron.canonize_info.is_flat = False
				new.Tetrahedron.canonize_info.face_status = {
					new.north_face(): tet.true_face_status(old.east_face()),
					new.south_face(): 2,
					new.east_face(): 2,
					new.west_face(): 2
					}
		if len(tet.Symmetries) == 12:
			new = self.new_arrow()
			old = Arrow(E02,F3,tet)
			new_copy = new.copy()
			for i in range(3):
				new.add_sym(new_copy)
				new_copy.rotate(1)
			# Face gluings.
			new.glue(new_copy.reverse())
			new.reverse().true_glue_as(old)
			# Edge labels.
			new.opposite()
			new.Tetrahedron.edge_labels = {
				new.equator(): tet.edge_labels[old.axis()],
				new.axis(): 3,
				new.north_head(): tet.edge_labels[old.north_head()],
				new.north_tail(): tet.edge_labels[old.south_head()],
				new.south_head(): 3,
				new.south_tail(): 3 
				}
			# Canonize info.
			if tet.canonize_info is not None:
				new.Tetrahedron.canonize_info = CanonizeInfo()
				new.Tetrahedron.canonize_info.part_of_coned_cell = True
				new.Tetrahedron.canonize_info.is_flat = False
				new.Tetrahedron.canonize_info.face_status = {
					new.north_face(): tet.true_face_status(old.east_face()),
					new.south_face(): 2,
					new.east_face(): 2,
					new.west_face(): 2
					}
		self.Tetrahedra.remove(tet)
		self.clear_and_rebuild()
		return


	def cancel_tetrahedra(self,edge):
		# First check valence and locus order.
		if edge.valence() == 2 and edge.LocusOrder == 1:
			case = 1
		elif edge.valence() == 1 and edge.LocusOrder == 2:
			case = 2
		else:
			return 0
		a = edge.get_arrow()
		if a.copy().next() is None:
			a.reverse()
		# In this non-geometric setting, we need to make sure the one-simplex opposite
		# edge has label 1 (ish). In the geometric setting, it actually might not be labelled
		# 1 and we can still do a cancellation. Need to be careful about this difference.
		# We will do this check within the two cases.
		#We don't want to try the cancellation if there are certain nontrivial symmetries.
		#They probably shouldn't be there anyway for a geometric triangulation, but we
		#check for them just in case.
		for sym in a.Tetrahedron.Symmetries:
			if sym.image(a.Edge) != a.Edge:
				return 0
		if case == 1:
			# valence 2, locus order 1.
			#Then the possible sub-cases are:
			#1. A single tet with faces properly glued to themselves, and without the symmetry.
			#2. Two tets without the symmetry.
			#3. Two tets both with the symmetry.
			if a.Tetrahedron.face_glued_to_self(a.Face):
				# Check the one-simplex opposite edge has label 1.
				if a.Tetrahedron.edge_labels[a.equator()] != 1:
					return 0
				#Check it's glued to itself properly.
				perm = a.Tetrahedron.Gluing[a.Face]
				if perm.image(a.Edge) != a.Edge:
					return 0
				b = a.copy().opposite()
				if len(a.Tetrahedron.Symmetries) != 1:
					# Don't think this can happen.
					return 0	
				c = b.copy().reverse()
				b.next().reverse()
				b.glue(c.glued())
				self.Tetrahedra.remove(a.Tetrahedron)
				self.clear_and_rebuild()
				new_edge = b.Tetrahedron.Class[b.Edge]	
				self.change_edge_labels(new_edge, 2)
			else:
				#In this case there are two distinct tetrahedra.
				b = a.glued()
				#Make sure the other tet doesn't have any illegal symmetries.
				for sym in b.Tetrahedron.Symmetries:
					if sym.image(b.Edge) != b.Edge:
						return 0
				#Another sanity check.
				if len(a.Tetrahedron.Symmetries) != len(b.Tetrahedron.Symmetries):
					return 0
				# Check that at least one of the two one-simplices opposite edge is labelled 1.
				if a.Tetrahedron.edge_labels[a.equator()] == 1: 
					new_label = b.Tetrahedron.edge_labels[b.equator()]
				elif b.Tetrahedron.edge_labels[b.equator()] == 1:
					new_label = a.Tetrahedron.edge_labels[a.equator()]
				else:
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
				self.clear_and_rebuild()
				new_edge = c.Tetrahedron.Class[c.Edge]	
				self.change_edge_labels(new_edge, new_label)
		if case == 2:
			# valence 1, locus order 2.
			# Then there are two sub-cases. 
			# 1. A single tet with faces adjacent to the edge glued to each other, with
			# a nontrivial symmetry taking the edge to itself. 
			# 2. A single tet with no non-trivial symmetries, the faces adjacent to the edge
			# glued to each other via a map which fixes the edge.
			if a.Tetrahedron.edge_labels[a.equator()] != 1:
				return 0
			b = a.copy().opposite()
			if len(a.Tetrahedron.Symmetries) == 2:
				# Sub-case 1.
				if b.copy().next() is None:
					b.reverse()
				b.next().reverse()
				b.glue(b.copy().reverse())
				self.Tetrahedra.remove(a.Tetrahedron)
				self.clear_and_rebuild()
				new_edge = b.Tetrahedron.Class[b.Edge]
				self.change_edge_labels(new_edge, 2)
			else:
				# Sub-case 2.
				b.next().reverse()
				b.glue(b.copy().reverse())
				b = a.copy().opposite().reverse()
				b.next().reverse()
				b.glue(b.copy().reverse())
				self.Tetrahedra.remove(a.Tetrahedron)
				self.clear_and_rebuild()
		return 1

	"""
	There are a few different cases of a 4-4 move. For now we just make the only one I think
	we need for canonize_part2. The version of this in CuspedOrbifold is called special_four_to_four.
	We use the same arrow diagrams as there.

	It takes three tets, one of which is flat with a sym, and makes it into 3 new tets, 2 with a sym.
	It's a kind of 4-4 move because you're re-triangulating an octahedron. This is very important
	in canonize_part2. It is used in step 2. Currently I don't have anything in this move to do with
	CanonizeInfo. I might add that at some point, or it might not be necessary.
	"""
	def four_to_four(self, edge):
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
		sym = a.Tetrahedron.nontrivial_sym()
		if sym.image(a.Edge) == a.Edge:
			return 0
		b = a.copy().true_next()
		c = b.glued()
		b_tet = b.Tetrahedron
		c_tet = c.Tetrahedron
		if b_tet is None or c_tet is None:
			return 0
		if b_tet is c_tet:
			return 0
		if len(b_tet.Symmetries) != 1 or len(c_tet.Symmetries) != 1:
			return 0
		# If we've gotten to here, then this case of a 4-4 move can be done.
		# We must determine if the symmetry axis is horizontal or vertical.
		if comp(a.south_face()) == sym.image(a.head()):
			case = 'horizontal'
		elif comp(a.north_face()) == sym.image(a.head()):
			case = 'vertical'
		else:
			raise Exception('error in simplicial four_to_four')
		new = self.new_arrows(3)
		if case == 'horizontal':
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
			new[0].Tetrahedron.edge_labels = {
				new[0].equator(): c_tet.edge_labels[c.south_head()],
				new[0].axis(): 1,
				new[0].north_head(): c_tet.edge_labels[c.equator()],
				new[0].north_tail(): c_tet.edge_labels[c.south_tail()],
				new[0].south_head(): c_tet.edge_labels[c.south_tail()],
				new[0].south_tail(): c_tet.edge_labels[c.equator()] 
				}
			new[1].Tetrahedron.edge_labels = {
				new[1].equator(): c_tet.edge_labels[c.north_head()],
				new[1].axis(): 1,
				new[1].north_head(): c_tet.edge_labels[c.north_tail()],
				new[1].north_tail(): c_tet.edge_labels[c.equator()],
				new[1].south_head(): b_tet.edge_labels[b.equator()],
				new[1].south_tail(): c_tet.edge_labels[c.south_tail()] 
				}
			new[2].Tetrahedron.edge_labels = {
				new[2].equator(): b_tet.edge_labels[b.north_tail()],
				new[2].axis(): 1,
				new[2].north_head(): b_tet.edge_labels[b.equator()],
				new[2].north_tail(): c_tet.edge_labels[c.north_tail()],
				new[2].south_head(): c_tet.edge_labels[c.north_tail()],
				new[2].south_tail(): b_tet.edge_labels[b.equator()] 
				}
		if case == 'vertical':
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

			new[0].Tetrahedron.edge_labels = {
				new[0].equator(): b_tet.edge_labels[b.south_tail()],
				new[0].axis(): 1,
				new[0].north_head(): c_tet.edge_labels[south_tail()],
				new[0].north_tail(): b_tet.edge_labels[b.equator()],
				new[0].south_head(): b_tet.edge_labels[b.equator()],
				new[0].south_tail(): c_tet.edge_labels[c.south_tail()] 
				}
			new[1].Tetrahedron.edge_labels = {
				new[1].equator(): c_tet.edge_labels[c.south_head()],
				new[1].axis(): 1,
				new[1].north_head(): c_tet.edge_labels[c.equator()],
				new[1].north_tail(): c_tet.edge_labels[c.south_tail()],
				new[1].south_head(): c_tet.edge_labels[c.north_tail()],
				new[1].south_tail(): b_tet.edge_labels[b.equator()] 
				}
			new[2].Tetrahedron.edge_labels = {
				new[2].equator(): c_tet.edge_labels[c.north_head()],
				new[2].axis(): 1,
				new[2].north_head(): c_tet.edge_labels[c.north_tail()],
				new[2].north_tail(): c_tet.edge_labels[c.equator()],
				new[2].south_head(): c_tet.edge_labels[c.equator()],
				new[2].south_tail(): c_tet.edge_labels[c.north_tail()] 
				}
		self.Tetrahedra.remove(a.Tetrahedron)
		self.Tetrahedra.remove(b.Tetrahedron)
		self.Tetrahedra.remove(c.Tetrahedron)
		self.clear_and_rebuild()
		return 1


	"""
	Sometimes we can collapse a tetrahedron onto one of its faces. We need this
	move in canonize_part2 when a face internal to the polyhedron is glued to itself.
	This can arise after a 1-4 move when one of the faces of the original tetrahedron
	is glued to itself.

	I don't currently have plans to use this anywhere other than canonize_part2. 
	"""
	def one_to_zero(self, two_subsimplex, tet):
		# Check that two_subsimplex is glued to itself.
		if tet.face_glued_to_self(two_subsimplex) is False:
			return 0
		# Now check that the "internal" edges have label 1. These are the edges of
		# tet which contain the vertex opposite two_subsimplex.
		inner_vertex = comp(two_subsimplex)
		for one_subsimplex in OneSubsimplices:
			if is_subset(inner_vertex,one_subsimplex) and tet.edge_labels[one_subsimplex] != 1:
				return 0
		# Find the edge which is mapped to itself by the face gluing map.
		perm = tet.Gluing[two_subsimplex]
		for one_subsimplex in OneSubsimplices:
			if (is_subset(one_subsimplex,two_subsimplex) and 
				perm.image(one_subsimplex) == one_subsimplex):
				break
		# We have two cases. When the symmetry group of tet has 3 elements,
		# the rotations of two_subsimplex, or when the symmetry group is trivial.
		# In either case, we start with this arrow.
		a = Arrow(one_subsimplex,two_subsimplex,tet)
		if len(tet.Symmetries) == 1:
			b = a.copy().opposite().reverse().next().reverse()
			c = a.copy().opposite().next()
			d = a.copy().reverse().next().reverse()
			# One final check: let's make sure that the tetrahedra these arrows belong
			# to are distinct.
			# (Actually this might not be necessary. Might change this.)
			tets_list = [e.Tetrahedron for e in (a,b,c,d)]
			if len(tets_list) != len(set(tets_list)):
				return 0
			b.glue(c)
			d.glue(d.copy().reverse())
			# Now we adjust the edge labels.
			b.Tetrahedron.edge_labels[b.axis()] = 2
			c.Tetrahedron.edge_labels[c.axis()] = 2
		elif len(tet.Symmetries) == 3 and tet.face_rotation(two_subsimplex):
			b = a.copy().reverse().true_next().reverse()
			if a.Tetrahedron is b.Tetrahedron:
				return 0
			b.glue(b.copy().reverse())
			b.Tetrahedron.edge_labels[b.north_head()] = 2
			b.Tetrahedron.edge_labels[b.south_head()] = 2
		else:
			return 0
		self.Tetrahedra.remove(tet)
		self.clear_and_rebuild()
		return 1

	"""
	The folowing method does the stellar move on a given edge. In other words, it introduces a
	finite vertex in the edge and cones the boundary of the star of the edge to that vertex.

	The tetrahedra around the edge are allowed to have the non-trivial order 2 symmetry which maps
	that edge to itself, but no other symmetries are allowed. We also require that any tet belonging
	to the star of the edge only has a single 1-simplex belonging to that edge's class (this actually
	implies the previous requirement about symmetries).
	"""
	def stellar_edge_move(self,edge):
		# Check if the move is possible.
		for corner in edge.Corners:
			tet = corner.Tetrahedron
			one_subsimplex = corner.Subsimplex
			for other_one_subsimplex in OneSubsimplices:
				if (other_one_subsimplex != one_subsimplex and 
					tet.Class[other_one_subsimplex] == tet.Class[one_subsimplex]):
					return 0
		# Now we know the move is possible.
		# Our next step it to find a tet with the symmetry or with a face
		# adjacent to edge which is glued to itself.
		a = edge.get_arrow()
		for i in range(len(edge.Corners)):
			if len(a.Tetrahedron.Symmetries) == 2:
				if a.glued().Tetrahedron is None:
					a.reverse()
				break
			if a.Tetrahedron.face_glued_to_self(a.Face):
				a.reverse()
				break
			a.reverse()
			if a.Tetrahedron.face_glued_to_self(a.Face):
				a.reverse()
				break
			a.reverse().next()
		# If a tet had the symmetry or a face glued to itself, then 'a' is now
		# in that tet. Otherwise, 'a' is just in some arbitrary tet containing
		# the edge.
		# Now we want to divide up a.Tetrahedron into two tets if it doesn't
		# have the symmetry, or one tet if it does.
		first_arrow = a.copy()
		if len(a.Tetrahedron.Symmetries) == 2:
			first_new = self.new_arrows(1)
			first_new[0].glue(first_new[0].copy().reverse())
			first_new[0].copy().reverse().true_glue_as(a.copy().opposite().reverse())
		else:
			first_new = self.new_arrows(2)
			first_new[0].glue(first_new[1].copy().reverse())
			first_new[0].copy().reverse().glue_as(a.copy().opposite().reverse())
			first_new[1].copy().reverse().glue_as(a.copy().opposite())
		extend_stellar_subdivision(a, first_new)
		return 1
		

	"""
	Helper function for the stellar_edge_move function above.

	The arrow a is in an old tet around the edge. new is a list with one or two elements,
	which are arrows for the new tets which are in a.Tetrahedron. We check if we've
	seen all tetrahedra already. If not, we continue it one more step into a.next(),
	and recurse.
	"""
	def extend_stellar_subdivision(self, a, new):
		tet = a.Tetrahedron
		# Assign the edge labels and canonize info for the new tets.
		b = a.copy()
		for c in new:
			c.Tetrahedron.edge_labels = {
				c.equator(): tet.edge_labels[b.axis()],
				c.axis(): tet.edge_labels[b.equator()],
				c.north_head(): 1,
				c.north_tail(): tet.edge_labels[b.south_head()],
				c.south_head(): 1,
				c.south_tail(): tet.edge_labels[b.south_tail()] 
				}
			if tet.canonize_info is not None:
				c.Tetrahedron.canonize_info = CanonizeInfo()
				c.Tetrahedron.canonize_info.part_of_coned_cell = True
				c.Tetrahedron.canonize_info.is_flat = False
				c.Tetrahedron.canonize_info.face_status = {
					c.north_face(): 2,
					c.south_face(): 2,
					c.east_face(): 2,
					c.west_face(): tet.canonize_info.face_status[b.south_face()]
					}
			b.reverse()
		# If either of the faces of tet adjacent to the edge are glued to themselves,
		# then we'll actually change some of the labels from 1 to 2. We do this later.
		# Now check if we've run through every tet in the star of the edge.
		# This will be true exactly when one of the following occurs.
		# 1. a.glued() == first_arrow.
		# 2. a.glued() == a.copy().reverse()
		# 3. a.glued().Tetrahedron is None
		# In (2), the face is glued to itself.
		# In (3), a.Tetrahedron has the symmetry and it is not the first tet.
		# If (2) occurs and a.Tetrahedron has the symmetry, then it must be
		# the first tet.
		if a.glued().Tetrahedron is None:
			self.Tetrahedra.remove(tet)
			self.clear_and_rebuild()
			return
		elif a.glued() == first_arrow:
			# In this case, it must be that there was no symmetry nor face glued to self.
			new[0].copy().opposite().glue(first_new[0].copy().opposite())
			new[1].copy().opposite().reverse().glue(first_new[1].copy().opposite().reverse())
			self.Tetrahedra.remove(tet)
			self.clear_and_rebuild()
			return
		elif a.glued() == a.copy().reverse():
			# Correct the edge labels.
			if len(tet.Symmetries) == 2:
				# Note this can only happen if a == first_arrow.
				# So the next case below is the other thing that can happen
				# if a == first_arrow.
				new[0].edge_labels[new[0].north_head()] = 2
				new[0].edge_labels[new[0].south_head()] = 2
				new[0].copy().opposite().glue(new[0].copy().opposite())
			elif a == first_arrow:
				# Then there is no sym, and both adjacent faces are glued to themselves.
				for c in new:
					c.edge_labels[c.north_head()] = 2
					c.edge_labels[c.south_head()] = 2
				new[0].copy().reverse().rotate(1).glue(new[1].copy().opposite().rotate(-1))
				new[1].copy().reverse().rotate(1).glue(new[0].copy().opposite().rotate(-1))
			else:
				new[0].edge_labels[north_head()] = 2
				new[1].edge_labels[south_head()] = 2
				new[0].copy().reverse().rotate(1).glue(new[1].copy().opposite().rotate(-1))
			self.Tetrahedra.remove(tet)
			self.clear_and_rebuild()
			return
		# If we've gotten to this point, then we know a.glued() lies in a new tet we must
		# extend to.
		b = a.glued()
		next_tet = b.Tetrahedron
		if len(next_tet.Symmetries) == 2:
			next_arrows = self.new_arrows(1)
			next_arrows[0].glue(next_arrows[0].copy().reverse())
			new[0].copy().opposite().glue(next_arrows[0].copy().opposite())
			next_arrows[0].copy().reverse().true_glue_as(b.copy().opposite().reverse())
			if len(tet.Symmetries) == 2:
				new[0].copy().opposite().reverse().glue(next_arrows[0].copy().opposite().reverse())
			else:
				new[1].copy().opposite().reverse().glue(next_arrows[0].copy().opposite().reverse())
		else:
			next_arrows = self.new_arrows(2)
			next_arrows[0].glue(next_arrows[1].copy().reverse())
			new[0].copy().opposite().glue(next_arrows[0].copy().opposite())
			if len(tet.Symmetries) == 2:
				new[0].copy().opposite().reverse().glue(next_arrows[1].copy().opposite().reverse())
			else:
				new[1].copy().opposite().reverse().glue(next_arrows[1].copy().opposite().reverse())
				next_arrows[0].copy().reverse().glue_as(b.copy().opposite().reverse())
				next_arrows[1].copy().reverse().glue_as(b.copy().opposite())
		self.extend_stellar_subdivision(b, next_arrows)
