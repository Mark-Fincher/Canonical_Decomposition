from tetrahedron import*
from arrow import*
from edge import*
from vertex import*
from simplex import*
from corner import*
from perm4 import*


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

	Also, we should allow for boundary faces now. It currently doesn't do that, but I'll change
	it. We should allow boundary faces because we might want to get the isomorphism signature
	of a simplicial orbifold with boundary, like if we're trying trying to build all triangulations
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
						for sym in corner.Tetrahedron:
							corner.Tetrahedron.Class[sym.image(corner.Subsimplex)] = newEdge
					#Now let's assign the LocusOrder of newEdge.
					newEdge.LocusOrder = corner.Tetrahedron.edge_labels[corner.Subsimplex]
		for i in range(len(self.Edges)):
			self.Edges[i].Index = i

	"""
	PACHNER MOVES.

	Here we define functions which do the 2-3, 3-2, 1-4, 4-1, 0-2, and 2-0 orbifold Pachner moves.
	They're similar to how they're written in the CuspedOrbifold class, except there's
	no geometry. And of course there are no 1-4 or 4-1 moves in the Cusped Orbifold class since
	they add or remove a non-ideal vertex.
	"""

	"""
	2-3 move. For now it does not check if a 2-3 move is valid to do. For instance, doesn't check if
	tet has the wrong kinds of symmetries or the face is glued to a different face of the same
	tet.
	"""
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
		# concisely, but I think it's easier to understand this way. 
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
			new.add_edge_label(3)
			new.opposite().add_edge_label(a.opposite().edge_label())
			new.opposite()
			a.opposite()
			for i in range(2):
				new.rotate(1).add_edge_label(a.rotate(1).edge_label())
			new.opposite().rotate(1)
			a.rotate(1)
			for i in range(2):
				new.rotate(1).add_edge_label(a.rotate(1).edge_label())
		elif tet.face_rotation(two_subsimplex):
			new = self.new_arrow()
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new.glue(new.copy())
			a.reverse()
			new.opposite()
			new.true_glue_as(a)
			new.reverse()
			new.true_glue_as(b)
			# Now we set the edge labels.
			new.add_edge_label(3)
			new.opposite().add_edge_label(b.opposite().edge_label())
			new.opposite()
			b.opposite()
			for i in range(2):
				new.rotate(1).add_edge_label(b.rotate(1).edge_label())
			new.rotate(1).reverse()
			for i in range(2):
				new.rotate(1).add_edge_label(a.rotate(1).edge_label())
			# Update canonize info.
			if tet.canonize_info is not None:
				new.Tetrahedron.canonize_info = CanonizeInfo()
				new.Tetrahedron.canonize_info.part_of_coned_cell = True
				new.Tetrahedron.canonize_info.is_flat = False
				new.Tetrahedron.canonize_info.face_status[new.Face] = 2
				new.Tetrahedron.canonize_info.face_status[new.copy().rotate(1).Face] = tet.canonize_info.face_status[a.copy().rotate(1).Face]
				new.Tetrahedron.canonize_info.face_status[new.copy().rotate(1).Face] = 2
				new.Tetrahedron.canonize_info.face_status[new.copy().rotate(1).reverse().Face] = b.Tetrahedron.canonize_info.face_status[b.copy().rotate(1).Face]
		elif tet.face_glued_to_self(two_subsimplex):
			new = self.new_arrows(2)
			for c in new:
				c.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			new[0].add_sym(new[0].copy().reverse())
			new[0].glue(new[1])
			new[1].glue(new[1].copy().reverse())
			# Now we glue the external faces.
			a.reverse()
			new[0].opposite()
			new[0].glue_as(a)
			a.rotate(-1)
			new[1].opposite()
			new[1].glue_as(a)
			a.rotate(-1)
			new[1].reverse()
			new[1].glue_as(a)
			# Now the edge labels.
			a.rotate(-1)
			new[0].add_edge_label(1)
			new[0].opposite().add_edge_label(a.opposite().edge_label())
			new[0].opposite()
			a.opposite()
			for i in range(2):
				new[0].rotate(1).add_edge_label(a.rotate(1).edge_label())
			new[0].opposite().rotate(1)
			a.rotate(1)
			for i in range(2):
				new[0].rotate(1).add_edge_label(a.rotate(1).edge_label())
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
				c.copy().opposite().add_edge_label(b.copy().opposite().edge_label())
				c.add_edge_label(1)
				c.rotate(-1).add_edge_label(b.rotate(-1).edge_label())
				c.rotate(-1).add_edge_label(b.rotate(-1).edge_label())
				c.rotate(1).reverse()
				c.rotate(1).add_edge_label(a.rotate(1).edge_label())
				c.rotate(1).add_edge_label(a.rotate(1).edge_label())
		self.Tetrahedra.remove(tet)
		if not tet.face_glued_to_self(two_subsimplex):
			# Then we need to remove b.Tetrahedron as well. If the face was glued to itself,
			# then there would be no b.
			self.Tetrahedra.remove(b.Tetrahedron)
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
			# See the diagram in the notes to make sense of the following.
			new = self.new_arrows(4)
			for i in range(4):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			for i in range(3):
				new[i].glue(new[(i+1)%3])
			new[3].glue(new[0].copy().opposite().rotate(-1))
			new[3].copy().rotate(-1).glue(new[2].copy().opposite().rotate(1))
			# Put an arrow in tet.
			a = Arrow(E12,F3,tet)
			new[0].copy().opposite().glue_as(a)
			new[3].copy().reverse().rotate(-1).glue_as(a.copy().reverse())
			new[1].copy().rotate(-1).glue_as(a.copy().reverse().rotate(1))
			new[2].copy().reverse().rotate(1).glue_as(a.copy().reverse().rotate(-1))
			# We're done with face gluings. Now edge labels. First deal with the central
			# edges, connected to the new vertex, which will all get labeled 1.
			central_arrows = [new[i].copy().opposite().reverse() for i in range(3)] + [new[3].copy().opposite()]
			for b in central_arrows:
				for i in range(3):
					b.rotate(1).add_edge_label(1)
			# Now the outer edges, with labels determined by the labels of tet.
			outer_new_arrows = [new[i].copy() for i in range(3)] + [new[3].copy().reverse()]
			outer_old_arrows = [a.copy().opposite(), a.copy().reverse().rotate(-1),
				a.copy().rotate(1).opposite(), a.copy().reverse().rotate(1)]
			# Now we just need to walk counter-clockwise around each face.
			for i in range(4):
				outer_new_arrows[i].add_edge_label(outer_old_arrows[i].edge_label())
				outer_old_arrows[i].rotate(1)
				outer_new_arrows[i].rotate(1)
				outer_new_arrows[i].add_edge_label(outer_old_arrows[i].edge_label())
				outer_old_arrows[i].reverse().rotate(1)
				outer_new_arrows[i].reverse().rotate(1)
				outer_new_arrows[i].add_edge_label(outer_old_arrows[i].edge_label())
			# Now update the canonize info.
			if tet.canonize_info is not None:
				old_arrows = [a.copy(),a.copy().rotate(-1),a.copy().rotate(1),a.copy().reverse()]
				new_out = [central.copy().reverse() for central in central_arrows]
				for i in range(4):
					new[i].Tetrahedron.canonize_info = CanonizeInfo()
					new[i].Tetrahedron.canonize_info.part_of_coned_cell = True
					new[i].Tetrahedron.canonize_info.is_flat = False
					for j in range(3):
						new[i].Tetrahedron.canonize_info.face_status[central_arrows[i].rotate(1).Face] = 2
					new[i].Tetrahedron.canonize_info.face_status[new_out[i].Face] = tet.canonize_info.face_status[old_arrows[i].Face] 
		if len(tet.Symmetries) == 2:
			new = self.new_arrows(2)
			for i in range(2):
				new[i].Symmetries.append(Perm4((0,1,2,3)))
			for sym in tet.Symmetries:
				if sym.tuple() != (0,1,2,3):
					break
			for one_subsimplex in OneSubsimplices:
				if sym.image(one_subsimplex) == one_subsimplex:
					a = Arrow(one_subsimplex,RightFace[one_subsimplex],tet)
					break
			a.opposite().rotate(-1)
			new[0].glue(new[1])
			new[0].copy().rotate(-1).glue(new[1].copy().reverse())
			b = new[0].copy().rotate(1).reverse()
			b.glue(b.copy().reverse())
			new[0].copy().opposite().reverse().true_glue_as(a)
			b = new[1].copy().opposite()
			b.glue(b.copy().reverse())
			b.reverse().true_glue_as(a.copy().rotate(1))
			b = new[0].copy().opposite()
			c = new[1].copy().opposite()
			for i in range(3):
				b.rotate(1).add_edge_label(1)
				c.rotate(1).add_edge_label(1)
			new[0].add_edge_label(a.copy().opposite())
			new[0].copy().rotate(-1).add_edge_label(a.copy().rotate(-1).edge_label())
			new[0].copy().reverse().rotate(1).add_edge_label(a.copy().rotate(1).edge_label())
			new[1].add_edge_label(a.copy().reverse().rotate(1).edge_label())
			new[1].copy().rotate(-1).add_edge_label(a.edge_label())
			new[1].copy().reverse().rotate(1).add_edge_label(a.copy().rotate(-1).edge_label())
			if tet.canonize_info is not None:
				for i in range(2):
					new[i].Tetrahedron.canonize_info = CanonizeInfo()
					new[i].Tetrahedron.canonize_info.part_of_coned_cell = True
					new[i].Tetrahedron.canonize_info.is_flat = False
				for i in range(3):
					new[0].Tetrahedron.canonize_info.face_status[b.rotate(1).Face] = 2
					new[1].Tetrahedron.canonize_info.face_status[c.rotate(1).Face] = 2
				new[0].Tetrahedron.canonize_info.face_status[new[0].copy().rotate(1).Face] = tet.canonize_info.face_status[a.Face]
				new[1].Tetrahedron.canonize_info.face_status[new[1].copy().rotate(1).Face] = tet.canonize_info.face_status[a.copy().rotate(1).Face]
		if len(tet.Symmetries) == 3:
			new = self.new_arrows(2)
			for i in range(2):
				new[i].Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			# Let Face be the face fixed by the rotations. Want to first get an
			# arrow with initial vertex the vertex fixed by rotations.
			for sym in tet.Symmetries:
				if sym.tuple() != (0,1,2,3):
					break
			for Face in TwoSubsimplices:
				if sym.image(Face) == Face:
					break
			for Edge in OneSubsimplices:
				if not is_subset(Edge,Face):
					a = Arrow(Edge,Face,tet)
					break
			new[0].glue(new[1])
			b = new[1].copy().reverse().rotate(1)
			b.glue(b)
			new[0].copy().opposite().reverse().glue_as(a)
			new[1].copy().reverse().rotate(-1).true_glue_as(a.copy().opposite().reverse())
			new[0].add_sym(new[0].copy().rotate(-1).reverse())
			new[0].add_sym(new[0].copy().reverse().rotate(1))
			# Now the edge labels.
			new[0].add_edge_label(a.copy().opposite().edge_label())
			new[0].copy().rotate(-1).add_edge_label(a.copy().rotate(-1).edge_label())
			new[0].copy().reverse().rotate(1).add_edge_label(a.copy().rotate(1).edge_label())
			b = new[0].copy().opposite()
			for i in range(3):
				b.rotate(1).add_edge_label(1)
			new[1].add_edge_label(a.copy().reverse().rotate(1))
			new[1].copy().rotate(-1).add_edge_label(a)
			new[1].copy().rotate(1).add_edge_label(a.copy().reverse().rotate(-1))
			b = new[1].copy().opposite()
			b.add_edge_label(1)
			b.rotate(1).add_edge_label(3)
			b.rotate(1).add_edge_label(1)
			if tet.canonize_info is not None:
				for i in range(2):
					new[i].Tetrahedron.canonize_info = CanonizeInfo()
					new[i].Tetrahedron.canonize_info.part_of_coned_cell = True
					new[i].Tetrahedron.canonize_info.is_flat = False
				b = new[0].copy().opposite()
				c = new[1].copy().opposite()
				for i in range(3):
					new[0].Tetrahedron.canonize_info.face_status[b.rotate(1).Face] = 2
					new[1].Tetrahedron.canonize_info.face_status[c.rotate(1).Face] = 2
				new[0].Tetrahedron.canonize_info.face_status[new[0].copy().rotate(1).Face] = tet.canonize_info.face_status[a.Face]
				new[1].Tetrahedron.canonize_info.face_status[new[1].copy().rotate(1).Face] = tet.canonize_info.face_status[a.copy().rotate(1).Face]
		if len(tet.Symmetries) == 4:
			new = self.new_arrow()
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			# Because of all the symmetries, we can take a to be any arrow in tet
			# and it will match our picture.
			a = Arrow(E12,F3,tet)
			b = new.copy().opposite().rotate(-1)
			b.glue(b.copy().reverse())
			b = new.copy().opposite()
			b.glue(b.copy().reverse())
			b = new.copy().rotate(1).reverse()
			b.glue(b.copy().reverse())
			new.copy().opposite().reverse().true_glue_as(a)
			b = new.copy().opposite()
			for i in range(3):
				b.rotate(1).add_edge_label(1)
			new.add_edge_label(a.copy().opposite())
			new.copy().rotate(-1).add_edge_label(a.copy().rotate(-1))
			new.copy().reverse().rotate(1).add_edge_label(a.copy().rotate(1))
			if tet.canonize_info is not None:
				new.Tetrahedron.canonize_info = CanonizeInfo()
				new.Tetrahedron.canonize_info.part_of_coned_cell = True
				new.Tetrahedron.canonize_info.is_flat = False
				for i in range(3):
					new.Tetrahedron.canonize_info.face_status[b.rotate(1).Face] = 2
				new.Tetrahedron.canonize_info.face_status[new.copy().rotate(1).Face] = tet.canonize_info.face_status[a.Face]
		if len(tet.Symmetries) == 12:
			new = self.new_arrow()
			new.Tetrahedron.Symmetries.append(Perm4((0,1,2,3)))
			a = Arrow(E12,F3,tet)
			b = new.copy().opposite().rotate(-1)
			b.glue(b.copy().reverse())
			new.copy().opposite().reverse().true_glue_as(a)
			b = new.copy().opposite()
			for i in range(3):
				b.rotate(1).add_edge_label(3)
			new.add_edge_label(a.copy().opposite())
			new.copy().rotate(-1).add_edge_label(a.copy().rotate(-1))
			new.copy().reverse().rotate(1).add_edge_label(a.copy().rotate(1))
			new.add_sym(new.copy().rotate(-1).reverse())
			new.add_sym(new.copy().reverse().rotate(1))
			if tet.canonize_info is not None:
				new.Tetrahedron.canonize_info = CanonizeInfo()
				new.Tetrahedron.canonize_info.part_of_coned_cell = True
				new.Tetrahedron.canonize_info.is_flat = False
				for i in range(3):
					new.Tetrahedron.canonize_info.face_status[b.rotate(1).Face] = 2
				new.Tetrahedron.canonize_info.face_status[new.copy().rotate(1).Face] = tet.canonize_info.face_status[a.Face]
		self.Tetrahedra.remove(tet)
		self.clear_and_rebuild()
		return


	def cancel_tetrahedra(self,edge):
		#First check valence and locus order.
		if edge.valence() == 2 and edge.LocusOrder == 1:
			case = 1
		elif edge.valence() == 1 and edge.LocusOrder == 2:
			case = 2
		else:
			return 0
		a = edge.get_arrow()
		if a.copy().next() is None:
			a.reverse()
		#After the cancellation, two edges will be collapsed on each other. We should
		#check that they have the same group label. They could actually be the same edge,
		#in which case they will of course have the same group label.
		if a.edge_label() != a.copy().next().edge_label():
			return 0
		#We don't want to try the cancellation if there are certain nontrivial symmetries.
		#They probably shouldn't be there anyway for a geometric triangulation, but we
		#check for them just in case.
		for sym in a.Tetrahedron.Symmetries:
			if sym.image(a.Edge) != a.Edge:
				return 0
		if case == 1:
			#Then the possible sub-cases are:
			#1. A single tet with faces properly glued to themselves, and without the symmetry.
			#2. Two tets without the symmetry.
			#3. Two tets both with the symmetry.
			#Split into two cases. If the faces adjacent to "edge" are glued to themselves,
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
				#In this case there are two distinct tetrahedra.
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
		if case == 2:
			#Then there's a single tet, with faces adjacent to the edge glued to each other. 
			#It could have the symmetry or not.
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
		self.clear_and_rebuild()
		return 1

	"""
	There are a few different cases of a 4-4 move. For now we just make the only one I think
	we need for canonize_part2.
	"""
	def four_to_four(self, edge):
		# This is the case where we have half of an octahedron, with the square face glued to 
		# itself by rotation around an axis which intersects the midpoints of two sides.
		if edge.valence() != 3 or edge.LocusOrder != 1:
			return 0
		a = edge.get_arrow()
		# First suppose a is in the tet with the symmetry (thought of as a flat tet).
		if len(a.Tetrahedron.Symmetries) == 2:
			if a.glued().Tetrahedron is None:
				a.opposite()
			# Now a is positioned where we want it.
		elif len(a.Tetrahedron.Symmetries) == 1:
			# Let's move the arrow into the tet with a nontrivial symmetry (if it exists).
			a.next()
			if len(a.Tetrahedron.Symmetries) != 2:
				a.next()
			if len(a.Tetrahedron.Symmetries) != 2:
				return 0
		else:
			return 0
		# Now the arrow a is positioned as in the diagram (in the write-up) and we do
		# the final checks.
		sym = a.Tetrahedron.nontrivial_sym()
		if sym.image(a.Edge) == a.Edge:
			return 0
		b = a.glued()
		c = b.glued()
		if b.Tetrahedron is None or c.Tetrahedron is None:
			return 0
		if b.Tetrahedron is c.Tetrahedron:
			return 0
		if len(b.Tetrahedron.Symmetries) != 1 or len(c.Tetrahedron.Symmetries) != 1:
			return 0
		# If we've gotten to here, then this case of a 4-4 move can be done.
		new = new_arrows(3)
		for i in range(3):
			new[i].Symmetries.append(Perm4((0,1,2,3)))
		for i in range(2):
			new[i].glue(new[i+1])
		new[0].add_sym(new[0].copy().reverse())
		new[2].add_sym(new[2].copy().reverse())
		new[0].copy().opposite().reverse().glue_as(b.copy().rotate(1))
		new[1].copy().opposite().glue_as(c.copy().reverse().rotate(-1))
		new[1].copy().opposite().reverse().glue_as(b.copy().rotate(-1))
		new[2].copy().opposite().reverse().glue_as(c.copy().reverse().rotate(1))
		# Now we add all the edge labels.
		for i in range(3):
			new[i].copy().opposite().add_edge_label(1)
		new[0].add_edge_label(b.copy().reverse().rotate(1).edge_label())
		new[0].copy().rotate(-1).add_edge_label(b.edge_label())
		new[0].copy().rotate(1).add_edge_label(b.copy().rotate(-1).edge_label())
		new[0].copy().reverse().rotate(1).add_edge_label(b.copy().rotate(-1).edge_label())
		new[0].copy().reverse().rotate(-1).add_edge_label(b.edge_label())
		new[1].add_edge_label(b.copy().reverse().rotate(-1).edge_label())
		new[1].copy().rotate(-1).add_edge_label(b.copy().rotate(1).edge_label())
		new[1].copy().rotate(1).add_edge_label(c.edge_label())
		new[1].copy().reverse().rotate(1).add_edge_label(b.edge_label())
		new[1].copy().reverse().rotate(-1).add_edge_label(b.copy().rotate(-1).edge_label())
		new[2].add_edge_label(c.copy().rotate(1).edge_label())
		new[2].copy().rotate(1).add_edge_label(b.copy().rotate(1).edge_label())
		new[2].copy().rotate(-1).add_edge_label(c.edge_label())
		new[2].copy().reverse().rotate(1).add_edge_label(b.copy().rotate(1).edge_label())
		new[2].copy().reverse().rotate(-1).add_edge_label(c.edge_label())
		self.Tetrahedra.remove(a.Tetrahedron)
		self.Tetrahedra.remove(b.Tetrahedron)
		self.Tetrahedra.remove(c.Tetrahedron)
		self.clear_and_rebuild()
		return 1



