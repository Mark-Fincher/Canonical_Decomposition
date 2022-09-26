"""
canonize_part2.py

For going from a proto-canonical triangulation to the canonical re-triangulation.

I'm following very closely the structure of the manifold case in SnapPy.
"""
from CanonizeInfo import*
from canonize_part1 import*
from SimplicialOrbifold import*
from IsomorphismSignature import hyperbolic_to_simplicial

"""
Input: orb, a CuspedOrbifold object which is proto-canonical, canonical decomp having cells
which are not tetrahedra.

Output: a SimplicialOrbifold object, the canonical retriangulation of the polyhedral canonical
decomposition of orb. It will have finite vertices, hence no longer geometric.
"""
def canonical_retriangulation(orb):
	# First, check if orb is already a canonical triangulation.
	if orb.is_proto_canonical() and not transparent_faces_or_flat_tets(orb):
		print("inputted orb is already a canonical triangulation")
		return orb
	# Check that it's proto-canonical.
	if not orb.is_proto_canonical():
		print("inputted orb is not proto-canonical")
		return
	# Set the canonize info.
	initialize_tet_status(orb)
	# Turn orb into a SimplicialOrbifold, i.e. remove the hyperbolic structure.
	# This function is in IsomorphismSignature.py.
	simplicial_orb = hyperbolic_to_simplicial(orb)
	step_one(simplicial_orb)
	step_two(simplicial_orb)
	# Assuming those two steps went okay, we now have the canonical retriangulation. Let's
	# remove the canonize_info data from the tetrahedra.
	for tet in simplicial_orb.Tetrahedra:
		tet.canonize_info = None
	return simplicial_orb

# Initializes the canonize info of each tet in the CuspedOrbifold orb. In particular,
# all the part_of_coned_cell flags are set to False.
def initialize_tet_status(orb):
	for tet in orb.Tetrahedra:
		new_canonize_info = CanonizeInfo()
		tet.canonize_info = new_canonize_info
		new_canonize_info.part_of_coned_cell = False
		new_canonize_info.is_flat = tet.is_flat()
		for face in TwoSubsimplices:
			if transparent_face(face,tet):
				new_canonize_info.face_status[face] = 1
			else:
				# Note: this will occur if tet or its neighbor along face is flat.
				# Or if face is glued to None.
				new_canonize_info.face_status[face] = 0

def step_one(orb):
	while cone_3_cell(orb) is True:
		pass

def cone_3_cell(orb):
	if insert_finite_vertex(orb) == 0:
		return False
	# Now expand the coned region.	
	while expand_coned_region(orb) == True:
		pass
	# If a transparent face is glued to itself, we do a 1-0 move on it.
	attempt_one_to_zero(orb)
	# attempt_cancellation is defined in canonize_part1.py.
	while attempt_cancellation(orb) == True:
		pass
	if verify_coned_region(orb) == False:
		raise Exception("in cone_3_cell, verify_coned_region returned False")
	return True

# NEW. To replace find_unconed_tet.
def insert_finite_vertex(orb):
	# First look for unconed tets with more than 2 syms. If there are any, do a 1-4 move on the tet
	# with the most syms.
	unconed_symmetric_tets = []
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell is False and tet.canonize_info.is_flat is False:
			if len(tet.Symmetries) > 2:
				unconed_symmetric_tets.append(tet)
	if len(unconed_symmetric_tets) > 0:
		# Find the tet with most syms from this list.
		tet_most_syms = unconed_symmetric_tets[0]
		for tet in unconed_symmetric_tets:
			if len(tet_most_syms.Symmetries) < len(tet.Symmetries):
				tet_most_syms = tet
		orb.one_to_four(tet_most_syms)
		return 1
	# 9/11/22 update. The following code replaces the commented out code beneath.
	# The point is that it's simpler than I was making it before:
	# the barycenter is in an edge if and only if it's transparent and
	# has locus order greater than 1.
	for edge in orb.Edges:
		if edge.LocusOrder > 1 and transparent_edge(edge):
			if orb.stellar_edge_move(edge):
				return 1
			else:
				raise Exception("tried and failed to do a stellar edge move")
	"""
	for edge in orb.Edges:
		if transparent_edge(edge) and more_than_one_rotation_intersecting(edge):
			if orb.stellar_edge_move(edge):
				return 1
			else:
				raise Exception("tried and failed to do a stellar edge move")
	"""
	# If we still haven't returned yet, then we now look for any unconed tet with
	# any non-trivial symmetry. That means it will have exactly 2 syms, since we
	# would have already coned it if it had more than 2.
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell is False and tet.canonize_info.is_flat is False:
			if len(tet.Symmetries) == 2:
				orb.one_to_four(tet)
				return 1
	# Finally, we'll do a 1-4 move on ANY unconed tet.
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell is False and tet.canonize_info.is_flat is False:
			orb.one_to_four(tet)
			return 1
	# If we haven't returned 1 by now, then there are no more 3-cells to cone.
	return 0

def transparent_edge(edge):
	# Returns True exactly when all faces adjacent to the edge are
	# transparent AND are not inside a coned polyhedron. The latter part
	# ensures that when we try to do a stellar edge move in the insert_finite_vertex
	# function, we only do it to an edge which is not already part of a coned poly.
	a = edge.get_arrow()
	for corner in edge.Corners:
		if a.Tetrahedron.true_face_status(a.Face) != 1:
			return False
		a.true_next()
	return True

# UPDATE 9/11/22. We don't need to use this function. The barycenter is in an edge if and
# only if the edge is transparent and has locus > 1. Checking if more than one rotation
# axis intersects is irrelevant. Will leave it here for now, but comment out the use of it,
# in insert_finite_vertex.
def more_than_one_rotation_intersecting(edge):
	# There are three ways a rotation axis could intersect an edge.
	# 1. The edge itself is a rotation axis, i.e. its locus order is > 1.
	# 2. A tet adjacent to the edge has the order 2 sym which is an involution of the edge.
	# 3. A face adjacent to the edge is glued to itself in such a way that the gluing
	#    map restricts to an involution of the edge.
	# We walk around the edge and count any occurences of these.
	a = edge.get_arrow()
	seen_tets = []
	if edge.LocusOrder > 1:
		num_axes = 1
	else:
		num_axes = 0
	for corner in edge.Corners:
		if a.Tetrahedron not in seen_tets:
			for sym in a.Tetrahedron.Symmetries:
				if sym.image(a.Edge) == a.Edge and sym.tuple() != (0,1,2,3):
					num_axes = num_axes + 1
					break
			b = a.copy().reverse()
			if a.glued() == b:
				num_axes = num_axes + 1
			if b.glued() == a:
				num_axes = num_axes + 1
				# We could have double counted an axis from this if the tet also
				# has the symmetry, but in that case there is certainly more than
				# 1 rotation axis intersecting the edge.
			if num_axes > 1:
				return True
		seen_tets.append(a.Tetrahedron)
		a.true_next()
	return False

def expand_coned_region(orb):
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell is True:
			for face in TwoSubsimplices:
				if tet.canonize_info.face_status[face] == 1:
					if tet.Neighbor[face].canonize_info.part_of_coned_cell is False:
						# Then we'll do a 2-3 move through face, as long as there are
						# either no symmetries or only rotations of face.
						if orb.check_two_to_three(face,tet) is False:
							if attempt_special_three_to_two(face,tet):
								return True
							else:
								raise Exception("unexpected symmetries of a tet while expanding cone region")
						if orb.two_to_three(face,tet):
							return True
						else:
							# This shouldn't happen.
							raise Exception("expand_coned_region, canonize_part2")
	return False


def attempt_special_three_to_two(face,tet):
	nbr = tet.Neighbor[face]
	if len(nbr.Symmetries) != 2:
		return 0
	for other_face in TwoSubsimplices:
		if tet.face_glued_to_self(other_face):
			one_subsimplex = other_face & face
			edge = tet.Class[one_subsimplex]
			if transparent_edge(edge) and edge.valence() == 3: 
				orb.three_to_two(edge)
				return 1
	return 0

def attempt_one_to_zero(orb):
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell:
			for face in TwoSubsimplices:
				if tet.canonize_info.face_status[face] == 1:
					if tet.face_glued_to_self(face):
						if orb.one_to_zero(face,tet):
							return
						else:
							raise Exception("error in attempt_one_to_zero")

def verify_coned_region(orb):
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell == True:
			for face in TwoSubsimplices:
				if tet.canonize_info.face_status[face] == 1:
					if tet.Neighbor[face].canonize_info.is_flat == False:
						return False
	return True

def step_two(orb):
	while eliminate_opaque_face(orb) == True:
		pass
	while attempt_cancellation(orb) == True:
		pass

def eliminate_opaque_face(orb):
	for tet in orb.Tetrahedra:
		# We only want to do a 2-3 or a 4-4 move to a tetrahedron which
		# is part of a coned cell. We don't want to do a move to one of the
		# diamonds this function has previously created which is not part of
		# a particular coned cell. However, a diamond will either not have
		# any canonize info, or all its faces will NOT have face_status opaque.
		# And in either case the tet will be skipped by this function.
		# We also want to skip flat tets, which exist just to express how
		# the face of a polyhedron is glued to itself.
		if tet.canonize_info is not None and tet.canonize_info.is_flat == False:
			for face in TwoSubsimplices:
				if tet.canonize_info.face_status[face] == 0:
					if tet.Neighbor[face].canonize_info.is_flat == False:
						if orb.check_two_to_three(face,tet) == False:
							raise Exception("error in eliminate_opaque_face")
						else:
							orb.two_to_three(face,tet)
							return True
					else:
						if opaque_four_to_four(face,tet) is False:
							raise Exception("error doing 4-4 through flat tet")
						return True
	return False

def opaque_four_to_four(face,tet):
	nbr = tet.Neighbor[face]
	if len(nbr.Symmetries) != 2:
		return False
	for one_subsimplex in OneSubsimplices:
		if is_subset(one_subsimplex,face):
			edge = tet.Class[one_subsimplex]
			if orb.four_to_four(edge):
				return True
	return False
				

