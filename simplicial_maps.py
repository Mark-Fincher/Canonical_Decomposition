# A simplicial map is just a map of orbifold triangulations. It sends tetrahedra
# to tetrahedra and respects the face gluings, symmetries, and edge labels. Hence
# it descends to a covering map of the quotient orbifold.

# Throughout this file, orb_1 and orb_2 can be geometric or not. Simplicial maps
# are combinatorial.

# Important convention: I represent a map from orb_1 to orb_2 as a dictionary

# {tet0: ((a, b, c, d), teti), ...}

# where tet0 of orb_1 is mapped to teti of orb_2 via the vertex mapping of vertices
# represented by (a,b,c,d), i.e. v0 -> va, v1 -> vb, v2 -> vc, v3 -> vd. Except
# it's not the tuple (a,b,c,d), but rather the perm4 object corresponding to it
# (see perm4.py). Similarly, tet0 and teti are the actual Tetrahedron objects,
# rather than just integers representing the index of the tet.

from SimplicialOrbifold import*
from HoroTriangle import* # for the function glued_to()

# See if orb_1 and orb_2 are simplicially isomorphic. This is true iff they cover each other.
def simplicial_isomorphic(orb_1,orb_2):
	if exists_covering(orb_1,orb_2) and exists_covering(orb_2,orb_1):
		return True
	return False

# See if there is a simplicial covering map from orb_1 to orb_2.
def exists_covering(orb_1,orb_2):
	if len(simplicial_maps(orb_1,orb_2)) > 0:
		return True
	return False

# See if orb_1 and orb_2 are OP simplicially isomorphic. This is true iff they OP cover each other.
def OP_simplicial_isomorphic(orb_1,orb_2):
	if exists_OP_covering(orb_1,orb_2) and exists_OP_covering(orb_2,orb_1):
		return True
	return False

# See if there is an OP simplicial covering map from orb_1 to orb_2.
def exists_OP_covering(orb_1,orb_2):
	if len(simplicial_maps_OP(orb_1,orb_2)) > 0:
		return True
	return False

# Find all simplicial maps from orb_1 to orb_2.
def simplicial_maps(orb_1,orb_2):
	simplicial_maps = []
	seen_maps_of_tet0 = []
	for tet in orb_2.Tetrahedra:
		for permutation in Perm4.S4():
			if (permutation.tuple(),tet) not in seen_maps_of_tet0:
				for sym in tet.Symmetries:
					seen_maps_of_tet0.append(((sym*permutation).tuple(),tet))
				simplicial_map = check_extends(orb_1,orb_2,permutation,tet)
				if simplicial_map != None:
					simplicial_maps.append(simplicial_map)
	return simplicial_maps

def simplicial_maps_OP(orb_1,orb_2):
	# Find the orientation preserving simplicial maps from orb_1 to orb_2.
	s_maps = simplicial_maps(orb_1,orb_2)
	s_maps_OP = [s_map for s_map in s_maps if s_map[orb_1.Tetrahedra[0]][0].sign() == 0]
	return s_maps_OP

"""
This checks if mapping tet_a to tet_b with map_a_to_b respects symmetries and 
the locus orders of the edges.
"""
def valid_tet_to_tet(tet_a,tet_b,map_a_to_b):
	for sym_a in tet_a.Symmetries:
		if (map_a_to_b*sym_a*inv(map_a_to_b)).tuple() not in [sym_b.tuple() for sym_b in tet_b.Symmetries]:
			return False
	for edge in OneSubsimplices:
		if tet_b.Class[map_a_to_b.image(edge)].LocusOrder % tet_a.Class[edge].LocusOrder != 0:
			return False
	return True

def check_extends(orb_1,orb_2,perm,image_of_tet0):
		simplicial_map = {orb_1.Tetrahedra[0]:(perm,image_of_tet0)}
		if valid_tet_to_tet(orb_1.Tetrahedra[0],image_of_tet0,perm) is False:
			return
		for tet in orb_1.Tetrahedra:
			if tet.Index > 0:
				# Create an entry in the mapping dictionary.
				simplicial_map[tet] = None
		active = [orb_1.Tetrahedra[0]]
		while active:
			tet = active.pop()
			image_of_tet = simplicial_map[tet][1]
			for face in TwoSubsimplices:
				voisin = tet.Neighbor[face]
				if voisin != None:
					# if we've already defined what the map should do to voisin, then we do nothing.
					# otherwise we now define what it does to voisin and add it to active.
					if simplicial_map[voisin] is None:
						phi = simplicial_map[tet][0]*inv(tet.Gluing[face])
						for sym in image_of_tet.Symmetries:
							if image_of_tet.Neighbor[sym.image(simplicial_map[tet][0].image(face))] != None:
								image_of_voisin = image_of_tet.Neighbor[sym.image(simplicial_map[tet][0].image(face))]
								simplicial_map[voisin] = (image_of_tet.Gluing[sym.image(simplicial_map[tet][0].image(face))]*sym*phi,image_of_voisin)
								if valid_tet_to_tet(voisin,image_of_voisin,simplicial_map[voisin][0]) is False:
									return
								break
						active.append(voisin)
		# check that simplicial_map respects the face gluings. By construction it respects some face gluings, but maybe not all.
		for tet1 in orb_1.Tetrahedra:
			for face1 in TwoSubsimplices:
				if tet1.Neighbor[face1] != None:
					tet2,face2 = glued_to(tet1,face1)
					well_defined_on_face1_and_face2 = False
					for sym1 in simplicial_map[tet1][1].Symmetries:
						for sym2 in simplicial_map[tet2][1].Symmetries:
							image_of_face1 = (sym1*simplicial_map[tet1][0]).image(face1)
							image_of_face2 = (sym2*simplicial_map[tet2][0]).image(face2)
							if simplicial_map[tet1][1].Neighbor[image_of_face1] != None:
								if glued_to(simplicial_map[tet1][1],image_of_face1) == (simplicial_map[tet2][1],image_of_face2):
									if (simplicial_map[tet1][1].Gluing[image_of_face1]*sym1*simplicial_map[tet1][0]).tuple() == (sym2*simplicial_map[tet2][0]*tet1.Gluing[face1]).tuple():
										well_defined_on_face1_and_face2 = True
					if well_defined_on_face1_and_face2 is False:
						return
		# if we haven't returned None by now, then simplicial_map is valid.
		return simplicial_map