"""
Take the quotient of a CuspedOrbifold under the action of a group of
simplicial automorphisms/self-isometries.

If the orbifold is geometric, the quotient will be too.

The function does not change the input orbifold. It returns a quotiented copy
of the input orbifold.

Each automorphism should be expressed as a dictionary of the form

{tet0: ((a, b, c, d), teti), ...}

The argument "group" should be a list or tuple of the full group you'd like to
quotient by, not just a set of generators of that group.
"""
from SimplicialOrbifold import*
from CuspedOrbifold import*

def quotient(initial_orb,initial_group):
	orb = initial_orb.copy()
	# initial_group refers to the tetrahedra of initial_orb.
	# We need to make a new list representing the group as it
	# acts on orb.
	group = []
	for initial_isom in initial_group:
		isom = {}
		for i in range(len(orb.Tetrahedra)):
			target_index = initial_isom[initial_orb.Tetrahedra[i]][1].Index
			isom[orb.Tetrahedra[i]] = (initial_isom[initial_orb.Tetrahedra[i]][0],
				orb.Tetrahedra[target_index])
		group.append(isom)
	orbits = make_orbits(orb,group)
	keep_tets = [orbit[0] for orbit in orbits]
	for tet in keep_tets:
		for face in TwoSubsimplices:
			update_face_gluing(tet,face,keep_tets,group)
	update_syms(orb,keep_tets,group)
	if not orb.is_geometric:
		# There could be some faces in a tet which are identified by a symmetry
		# but more than one of them is not glued to None. At this point, that can only
		# happen if each is glued to itself. We don't allow that for simplicial orbs.
		for tet in keep_tets:
			for face in TwoSubsimplices:
				if tet.face_glued_to_self(face):
					for sym in tet.Symmetries:
						other_face = sym.image(face)
						if other_face != face:
							tet.Neighbor[other_face] = None
							tet.Gluing[other_face] = None
		# Now we have to be careful about edge labels.
		# In the geometric case, this is taken care of by the geometry.
		# The point is that if an edge has a non-trivial stabilizer under the
		# action of group, then its new edge label will be the old one times the
		# order of that stabilizer. The order of the stabilizer is the valence of
		# the edge in the original orbifold divided by the valence in the quotient
		# orbifold. So we want to remember the old valences before we destroy them
		# in making the new edge classes.
		old_edge_class = {}
		for tet in keep_tets:
			for one_subsimplex in OneSubsimplices:
				old_edge_class[(tet,one_subsimplex)] = tet.Class[one_subsimplex]
			tet.clear_Class()
		quotient_orb = SimplicialOrbifold(keep_tets)
		seen_new_edges = []
		for tet in quotient_orb.Tetrahedra:
			for one_subsimplex in OneSubsimplices:
				new_edge = tet.Class[one_subsimplex]
				if new_edge in seen_new_edges:
					continue
				seen_new_edges.append(new_edge)
				old_edge = old_edge_class[(tet,one_subsimplex)]
				stabilizer_order = old_edge.valence()//new_edge.valence()
				new_label = old_edge.LocusOrder*stabilizer_order
				quotient_orb.change_edge_labels(new_edge,new_label)
		return quotient_orb
	# When there's geometry, the edge labels will be determined by that geometry,
	# so we don't need to do the same work as above. Might still be good to do it,
	# as a check.
	for tet in keep_tets:
		tet.remove_redundant_face_gluing()
		tet.clear_Class()
	
	for i in range(len(keep_tets)):
			print('gluing data for', keep_tets[i], 'is')
			for j in range(4):
				print('face', j, 'of', keep_tets[i], 'glued to', keep_tets[i].Neighbor[TwoSubsimplices[j]])
				print('with gluing map',keep_tets[i].Gluing[TwoSubsimplices[j]])
			print('symmetries of',keep_tets[i],'are')
			for sym in keep_tets[i].Symmetries:
				print(sym)

	return CuspedOrbifold(keep_tets)

def quotient_old(orb,group):
	orbits = make_orbits(orb,group)
	keep_tets = []
	for orbit in orbits:
		keep_tets.append(orbit[0])
	for tet in keep_tets:
		for face in TwoSubsimplices:
			update_face_gluing(tet,face,keep_tets,group)
	update_syms(orb,keep_tets,group)
	if type(orb) is SimplicialOrbifold:
		# There could be some faces in a tet which are identified by a symmetry
		# but more than one of them is not glued to None. At this point, that can only
		# happen if each is glued to itself. We don't allow that for simplicial orbs.
		for tet in keep_tets:
			for face in TwoSubsimplices:
				if tet.face_glued_to_self(face):
					for sym in tet.Symmetries:
						other_face = sym.image(face)
						if other_face != face:
							tet.Neighbor[other_face] = None
							tet.Gluing[other_face] = None
		# Now we have to be careful about edge labels.
		# In the CuspedOrbifold case, this is taken care of by the geometry.
		# The point is that if an edge has a non-trivial stabilizer under the
		# action of group, then its new edge label will be the old one times the
		# order of that stabilizer. The order of the stabilizer is the valence of
		# the edge in the original orbifold divided by the valence in the quotient
		# orbifold. So we want to remember the old valences before we destroy them
		# in making the new edge classes.
		old_edge_class = {}
		for tet in keep_tets:
			for one_subsimplex in OneSubsimplices:
				old_edge_class[(tet,one_subsimplex)] = tet.Class[one_subsimplex]
			tet.clear_Class()
		quotient_orb = SimplicialOrbifold(keep_tets)
		seen_new_edges = []
		for tet in quotient_orb.Tetrahedra:
			for one_subsimplex in OneSubsimplices:
				new_edge = tet.Class[one_subsimplex]
				if new_edge in seen_new_edges:
					continue
				seen_new_edges.append(new_edge)
				old_edge = old_edge_class[(tet,one_subsimplex)]
				stabilizer_order = old_edge.valence()//new_edge.valence()
				new_label = old_edge.LocusOrder*stabilizer_order
				quotient_orb.change_edge_labels(new_edge,new_label)
		return quotient_orb
	# The CuspedOrbifold case is simpler.The only thing we need to possibly
	# adjust at this point is that some faces might be glued to None when they
	# could be glued to themselves. This is against my convention for CuspedOrbifolds,
	# which is of course different from my convention for SimplicialOrbifolds. Eventually
	# I will have the same conventions for each. Dealing with other things now.
	for tet in keep_tets:
		for face in TwoSubsimplices:
			if tet.Neighbor[face] is None:
				for sym in tet.Symmetries:
					other_face = sym.image(face)
					if tet.face_glued_to_self(other_face):
						gluing = tet.Gluing[other_face]
						tet.attach(face,tet,(inv(sym)*gluing*sym).tuple())
						break
		tet.clear_Class()
	return CuspedOrbifold(keep_tets)

def make_orbits(orb,group):
	orbits_list = []
	seen_tets = []
	for tet in orb.Tetrahedra:
		if tet not in seen_tets:
			seen_tets.append(tet)
			new_orbit = [tet]
			for g in group:
				image_tet = g[tet][1]
				if image_tet not in seen_tets:
					seen_tets.append(image_tet)
				if image_tet not in new_orbit:
					new_orbit.append(g[tet][1])
			orbits_list.append(new_orbit)
	return orbits_list

def update_face_gluing(tet,face,keep_tets,group):
	nbr = tet.Neighbor[face]
	# If nbr is None, do nothing.
	if nbr is None:
		return
	if nbr not in keep_tets:
		# Find the image of nbr which IS in keep_tets.
		for g in group:
			if g[nbr][1] in keep_tets:
				break
		gluing = tet.Gluing[face]
		tet.attach(face,g[nbr][1],(g[nbr][0]*gluing).tuple())
		# Now, if there is some g in the group such that g(face) != face, we want to
		# remove the face gluing data for g(face). That's just our convention for
		# orbifold triangulations. But it's also useful here because
		# it means update_face_gluing will skip g(face) when it gets to it later.
		for g in group:
			if g[tet][1] is tet:
				other_face = g[tet][0].image(face)
				if other_face != face:
					tet.Neighbor[other_face] = None
					tet.Gluing[other_face] = None
	return

# Note: to update sym group for each tet, we don't just add the elements of the isom group
# which fix the tet, rather the group they generate with the original sym group of the tet. (?)
def update_syms(orb,keep_tets,group):
	for tet in keep_tets:
		generators = tet.Symmetries
		for g in group:
			if g[tet][1] == tet and g[tet][0].tuple() != (0,1,2,3):
				generators.append(g[tet][0])
		tet.Symmetries = Perm4.GenerateSubgroup(generators)
	# If F and F' are faces in a tet which are identified by a face gluing
	# and by a symmetry, then they should really be glued to themselves.
	# This can't happen in the original orbifold, but it could happen now
	# that we've possibly added more symmetries. Update: this face gluing convention
	# isn't something I believe in anymore, should change this at some point.
	for tet in keep_tets:
		for face1 in TwoSubsimplices:
			if tet.Neighbor[face1] is tet:
				gluing = tet.Gluing[face1]
				face2 = gluing.image(face1)
				if face1 != face2:
					for sym in tet.Symmetries:
						if sym.image(face1) == face2:
							tet.attach(face1,tet,(inv(gluing)*sym).tuple())
							tet.attach(face2,tet,(gluing*inv(sym)).tuple())
							break
