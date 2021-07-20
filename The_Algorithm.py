"""
This file contains the canonical decomposition algorithm. Each function in the algorithm is
a method in the CuspedOrbifold class (or not?).
"""
from Exact_Arithmetic import*
from simplex import*
from CuspedOrbifold import*
from edge import*
from HoroTriangle import*
from tetrahedron import*
from vertex import*
from Dest_to_Triang import*
#from sagestuff import*

"""
Need interval arithmetic for this.
"""

def check_2_to_3_possible(tets,tet,face):
	for edge in OneSubsimplices:
		if is_subset(edge,face):
			z = tet.edge_params[edge]
			w = tet.Neighbor[face].edge_params[tet.Gluing[face].image(edge)]
			if (z*w).imag.evaluate() < 0 or (z*w).imag = SquareRootCombination.Zero():
				return False
	# Now we check that tet doesn't have symmetries taking "face" to a a different face. If so, then a 2-3 move won't make
	# sense in the universal cover.
	for perm in tet.Symmetries:
		if perm[FaceIndex[face]] != FaceIndex[face]:
			return False
	for perm in tet.Neighbor[face].Symmetries:
		if perm[tet.Gluing[face][FaceIndex[face]]] != tet.Gluing[face][FaceIndex[face]]:
			return False
	# If we haven't returned False by now, then return True.
	return True



def two_to_three(triang,tet,face):
	#currently triang should be a list of tets, not a CuspedOrbifold object. That could change.
	#First we determine if the face is glued to itself and/or tet has rotational symmetries fixing
	#the vertex FaceIndex[face]. To understand this function, assume they're both false at first
	#and read the write-up.
	if tet.Neighbor[face] == tet:
		face_glued_to_self = True
		for i in range(4):
			if tet.Gluing[face][i] == i:
				v0 = i
	else:
		face_glued_to_self = False
	fix_face_perms = []
	for perm in tet.Symmetries:
		if perm[FaceIndex[face]] == FaceIndex[face]:
			fix_face_perms.append(perm)
	if len(fix_face_perms) > 1:
		face_rotation = True
	else:
		face_rotation = False
	custom_PickAnEdge = {F0:(3,1),F1:(2,0),F2:(0,1),F3:(2,1)}
	a1,b1 = custom_PickAnEdge[face]
	d1 = FaceIndex[face]
	c1 = ({0,1,2,3}-{a1,b1,d1}).pop()
	z0 = tet.edge_params[bitmap((a1,b1))]
	if face_glued_to_self:
		# Then we make a new tetrahedron which is a copy of tet, glue it to tet along face		
		perm = tet.Gluing[face]
		tet.detach(face)
		other_tet = Tetrahedron()
		other_tet.edge_params = tet.edge_params
		other_tet.Symmetries = tet.Symmetries
		tet.attach(face,other_tet,[perm[0],perm[1],perm[2],perm[3]])
	else:
		other_tet = tet.Neighbor[face]
	a2,b2,c2,d2 = tet.Gluing[face]((a1,b1,c1,d1))
	w0 = other_tet.edge_params[bitmap((c2,b2))]
	new_tet0 = Tetrahedron()
	new_tet1 = Tetrahedron()
	new_tet2 = Tetrahedron()
	new_tet0.attach(F0,new_tet1,[1,0,2,3])
	if tet.Neighbor[TwoSubsimplices[b1]] != None:
		perm = tet.Gluing[TwoSubsimplices[b1]]
		voisin = tet.Neighbor[TwoSubsimplices[b1]]
		tet.detach(TwoSubsimplices[b1])
		new_tet0.attach(F1,voisin,[perm[a1],perm[b1],perm[d1],perm[c1]])
	if other_tet.Neighbor[TwoSubsimplices[b2]] != None:
		perm = other_tet.Gluing[TwoSubsimplices[b2]]
		voisin = other_tet.Neighbor[TwoSubsimplices[b2]]
		other_tet.detach(TwoSubsimplices[b2])
		new_tet0.attach(F2,voisin,[perm[a2],perm[d2],perm[b2],perm[c2]])
	new_tet0.attach(F3,new_tet2,[2,0,3,1])
	if tet.Neighbor[TwoSubsimplices[a1]] != None:
		perm = tet.Gluing[TwoSubsimplices[a1]]
		voisin = tet.Neighbor[TwoSubsimplices[a1]]
		tet.detach(TwoSubsimplices[a1])
		new_tet1.attach(F0,voisin,[perm[a1],perm[b1],perm[d1],perm[c1]])
	if other_tet.Neighbor[TwoSubsimplices[a2]] != None:
		perm = other_tet.Gluing[TwoSubsimplices[a2]]
		voisin = other_tet.Neighbor[TwoSubsimplices[a2]]
		other_tet.detach(TwoSubsimplices[a2])
		new_tet1.attach(F2,voisin,[perm[d2],perm[b2],perm[a2],perm[c2]])
	new_tet1.attach(F3,new_tet2,[0,1,3,2])
	if tet.Neighbor[TwoSubsimplices[c1]] != None:
		perm = tet.Gluing[TwoSubsimplices[c1]]
		voisin = tet.Neighbor[TwoSubsimplices[c1]]
		tet.detach(TwoSubsimplices[c1])
		new_tet2.attach(F0,voisin,[perm[c1],perm[b1],perm[a1],perm[d1]])
	if other_tet.Neighbor[TwoSubsimplices[c2]] != None:
		perm = other_tet.Gluing[TwoSubsimplices[c2]]
		voisin = other_tet.Neighbor[TwoSubsimplices[c2]]
		other_tet.detach(TwoSubsimplices[c2])
		new_tet2.attach(F3,voisin,[perm[d2],perm[b2],perm[a2],perm[c2]])
	new_tet0.fill_edge_params(z0*w0/(z0 - ComplexSquareRootCombination.One() + w0))
	new_tet1.fill_edge_params((z0 - ComplexSquareRootCombination.One() + w0)/w0)
	new_tet2.fill_edge_params((w0 - ComplexSquareRootCombination.One())/(z0 - ComplexSquareRootCombination.One() + w0))
	if face_glued_to_self and not face_rotation:
		# v0 is the vertex fixed by the face gluing map
		if v0 == b1:	
			new_tet0.Symmetries = [Perm4((0,1,2,3)),Perm4((3,2,1,0))]
			new_tet0.detach(F3)
			new_tet1.detach(F3)
			new_tet1.attach(F3,new_tet1,[2,1,0,3])
			new_tets = [new_tet0,new_tet1]
			for T in triang:
				if T != tet:
					new_tets.append(T)
			index_tets(new_tets)
			return new_tets
		elif v0 == a1:
			new_tet1.Symmetries = [Perm4((0,1,2,3)),Perm4((2,3,0,1))]
			new_tet1.detach(F3)
			new_tet0.detach(F3)
			new_tet0.attach(F3,new_tet0,[0,2,1,3])
			new_tets = [new_tet0,new_tet1]
			for T in triang:
				if T != tet:
					new_tets.append(T)
			index_tets(new_tets)
			return new_tets
		elif v0 == c1:
			new_tet2.Symmetries = [Perm4((0,1,2,3)),Perm4((3,2,1,0))]
			new_tet2.detach(F2)
			new_tet0.detach(F0)
			new_tet0.attach(F0,new_tet0,[0,2,1,3])
			new_tets = [new_tet0,new_tet2]
			for T in triang:
				if T != tet:
					new_tets.append(T)
			index_tets(new_tets)
			return new_tets
	elif face_rotation and not face_glued_to_self:
		new_tet0.detach(F3)
		new_tet0.detach(F0)
		new_tet0.attach(F0,new_tet0,[3,1,2,0])
		new_tets = [new_tet0]
		for T in triang:
			if T != tet and T != other_tet:
				new_tets.append(T)
		index_tets(new_tets)
		return new_tets
	elif face_glued_to_self and face_rotation:
		new_tet0.detach(F3)
		new_tet0.detach(F0)
		new_tet0.attach(F0,new_tet0,[3,1,2,0])
		new_tet0.Symmetries = [Perm4((0,1,2,3)),Perm4((3,2,1,0))]
		new_tets = [new_tet0]
		for T in triang:
			if T != tet:
				new_tets.append(T)
		index_tets(new_tets)
		return new_tets
	else:
		new_tets = [new_tet0,new_tet1,new_tet2]
		for T in triang:
			if T != tet and T != other_tet:
				new_tets.append(T)
		index_tets(new_tets)
		return new_tets





