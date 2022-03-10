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
	# Now update the CanonizeInfo.
	for tet in orb.Tetrahedra:
		new_canonize_info = CanonizeInfo()
		tet.canonize_info = new_canonize_info
		new_canonize_info.part_of_coned_cell = False
		new_canonize_info.is_flat = tet.is_flat()
		for face in TwoSubsimplices:
			if transparent_face(face,tet):
				new_canonize_info.face_status[face] = 1
			else:
				new_canonize_info.face_status[face] = 0
	# Turn orb into a SimplicialOrbifold.
	orb = hyperbolic_to_simplicial(orb)
	step_one(orb)
	step_two(orb)
	# Assuming those two steps went okay, we now have the canonical retriangulation. Let's
	# remove the canonize_info data from the tetrahedra.
	for tet in orb.Tetrahedra:
		tet.canonize_info = None
	return 1

def step_one(orb):
	num_finite_vertices = 0
	while cone_3_cell(orb,num_finite_vertices) is True:
		pass

def cone_3_cell(orb,num_finite_vertices):
	# Find an unconed_tet with the most symmetries among all such tets.
	unconed_tets = []
	for tet in orb.Tetrahedra:
		if tet.canonize_info.part_of_coned_cell is False and tet.canonize_info.is_flat is False:
			unconed_tets.append(tet)
	if len(unconed_tets) == 0:
		return False
	tet_most_syms = unconed_tets[0]
	for tet in unconed_tets:
		if len(tet_most_syms.Symmetries) < len(tet.Symmetries):
			tet_most_syms = tet
	orb.one_to_four(tet_most_syms)

def step_two(orb):
	# more stuff


