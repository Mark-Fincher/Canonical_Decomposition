"""
canonize_part2.py

For going from a proto-canonical triangulation to the canonical re-triangulation.

I'm following very closely the structure of the manifold case in SnapPy.
"""
from CanonizeInfo import*
from canonize import*
from SimplicialOrbifold import*
from IsomorphismSignature import hyperbolic_to_simplicial

"""
Input: orb, a CuspedOrbifold object which is proto-canonical, canonical decomp having cells
which are not tetrahedra.

Output: a SimplicialOrbifold object, the canonical retriangulation or the polyhedral canonical
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
				new_canonize_info.face_status[face] = 0

def step_one(orb):
	while cone_3_cell(orb) is True:
		print('cone 3 cell was true')
		pass

def cone_3_cell(orb):
	tet = find_unconed_tet(orb)
	if tet is None:
		return False
	"""
	orb.one_to_four(tet)
	# Now we've introduced a finite vertex. Want it to be the barycenter of the polyhedron.
	# We might need to collapse some of the newly created tetrahedra so that this finite
	# vertex really is at the barycenter.
	move_to_barycenter(orb)
	"""
	# Instead of the above, maybe I do some kind of insert_finite_vertex(tet).
	# Now expand the coned region.	
	while expand_coned_region(orb) == True:
		print('expanded cone region')
		pass
	# attempt_cancellation is defined in canonize.py.
	while attempt_cancellation(orb) == True:
		print('cancelled tet(s)')
		pass
	if verify_coned_region(orb) == False:
		raise Exception("in cone_3_cell, verify_coned_region returned False")
	return True

# Find an unconed tet. Prioritize tets with transparent faces which are glued to themselves or with
# the max number of symmetries over all choices. If all tets are coned, return None.
def find_unconed_tet(orb):
	unconed_tets = []
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell is False and tet.canonize_info.is_flat is False:
			for face in TwoSubsimplices:
				if tet.canonize_info.face_status[face] == 1 and tet.face_glued_to_self(face):
					return tet
			unconed_tets.append(tet)
	if len(unconed_tets) == 0:
		return
	tet_most_syms = unconed_tets[0]
	for tet in unconed_tets:
		if len(tet_most_syms.Symmetries) < len(tet.Symmetries):
			tet_most_syms = tet
	return tet_most_syms


# Probably replacing this with insert_finite_vertex somehow
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

def expand_coned_region(orb):
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell is True:
			for face in TwoSubsimplices:
				if tet.canonize_info.face_status[face] == 1:
					if tet.Neighbor[face].canonize_info.part_of_coned_cell is False:
						# Then we'll do a 2-3 move through face, as long as there are
						# either no symmetries or only rotations of face. And if face
						# is not glued to itself.
						if tet.face_glued_to_self(face):
							raise Exception("unexpected face glued to self while expanding cone region")
						if orb.check_two_to_three(face,tet) is False:
							raise Exception("unexpected symmetries of a tet while expanding cone region")
						if orb.two_to_three(face,tet):
							return True
						else:
							# This shouldn't happen.
							raise Exception("expand_coned_region, canonize_part2")
	return False

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
				

