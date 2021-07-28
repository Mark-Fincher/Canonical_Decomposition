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
from sagestuff import*

"""
Need interval arithmetic for this.
"""

def check_2_to_3_possible(tets,tet,face):
	for edge in OneSubsimplices:
		if is_subset(edge,face):
			z = tet.edge_params[edge]
			w = tet.Neighbor[face].edge_params[tet.Gluing[face].image(edge)]
			if (z*w).imag.evaluate() < 0 or (z*w).imag == SquareRootCombination.Zero():
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
	if tet.Neighbor[face] == tet and tet.Gluing[face][FaceIndex[face]] == FaceIndex[face]:
		face_glued_to_self = True
		for i in range(4):
			if tet.Gluing[face][i] == i and i != FaceIndex[face]:
				vert0 = i
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
	# Before we get started assigning all the face pairings, we define two useful dictionaries. See the write-up
	# for what they're about.
	perm_key = {(tet,b1):Perm4((a1,b1,d1,c1)),(tet,a1):Perm4((a1,b1,d1,c1)),(tet,c1):Perm4((c1,b1,a1,d1)),(other_tet,b2):Perm4((a2,d2,b2,c2)),(other_tet,a2):Perm4((d2,b2,a2,c2)),(other_tet,c2):Perm4((d2,b2,a2,c2))}
	tet_key = {(tet,b1):new_tet0,(tet,a1):new_tet1,(tet,c1):new_tet2,(other_tet,b2):new_tet0,(other_tet,a2):new_tet1,(other_tet,c2):new_tet2} 
	# Now we do the face pairings.
	new_tet0.attach(F0,new_tet1,[1,0,2,3])
	if tet.Neighbor[TwoSubsimplices[b1]] != None:
		perm = tet.Gluing[TwoSubsimplices[b1]]
		voisin = tet.Neighbor[TwoSubsimplices[b1]]
		key1 = (tet,b1)
		if voisin == tet or voisin == other_tet:
			key2 = (voisin,FaceIndex[perm.image(TwoSubsimplices[b1])])
			new_tet0.attach(F1,tet_key[key2],(inv(perm_key[key2])*perm*perm_key[key1]).tuple())
		else:
			tet.detach(TwoSubsimplices[b1])
			new_tet0.attach(F1,voisin,(perm*perm_key[key1]).tuple())
	if other_tet.Neighbor[TwoSubsimplices[b2]] != None:
		perm = other_tet.Gluing[TwoSubsimplices[b2]]
		voisin = other_tet.Neighbor[TwoSubsimplices[b2]]
		key1 = (other_tet,b2)
		if voisin == tet or voisin == other_tet:
			key2 = (voisin,FaceIndex[perm.image(TwoSubsimplices[b2])])
			new_tet0.attach(F2,tet_key[key2],(inv(perm_key[key2])*perm*perm_key[key1]).tuple())
		else:
			other_tet.detach(TwoSubsimplices[b2])
			new_tet0.attach(F2,voisin,(perm*perm_key[key1]).tuple())
	new_tet0.attach(F3,new_tet2,[2,0,3,1])
	if tet.Neighbor[TwoSubsimplices[a1]] != None:
		perm = tet.Gluing[TwoSubsimplices[a1]]
		voisin = tet.Neighbor[TwoSubsimplices[a1]]
		key1 = (tet,a1)
		if voisin == tet or voisin == other_tet:
			key2 = (voisin,FaceIndex[perm.image(TwoSubsimplices[a1])])
			new_tet1.attach(F0,tet_key[key2],(inv(perm_key[key2])*perm*perm_key[key1]).tuple())
		else:
			tet.detach(TwoSubsimplices[a1])
			new_tet1.attach(F0,voisin,(perm*perm_key[key1]).tuple())
	if other_tet.Neighbor[TwoSubsimplices[a2]] != None:
		perm = other_tet.Gluing[TwoSubsimplices[a2]]
		voisin = other_tet.Neighbor[TwoSubsimplices[a2]]
		key1 = (other_tet,a2)
		if voisin == tet or voisin == other_tet:
			key2 = (voisin,FaceIndex[perm.image(TwoSubsimplices[a2])])
			new_tet1.attach(F2,tet_key[key2],(inv(perm_key[key2])*perm*perm_key[key1]).tuple())
		else:
			other_tet.detach(TwoSubsimplices[a2])
			new_tet1.attach(F2,voisin,(perm*perm_key[key1]).tuple())
	new_tet1.attach(F3,new_tet2,[0,1,3,2])
	if tet.Neighbor[TwoSubsimplices[c1]] != None:
		perm = tet.Gluing[TwoSubsimplices[c1]]
		voisin = tet.Neighbor[TwoSubsimplices[c1]]
		key1 = (tet,c1)
		if voisin == tet or voisin == other_tet:
			key2 = (voisin,FaceIndex[perm.image(TwoSubsimplices[c1])])
			new_tet2.attach(F0,tet_key[key2],(inv(perm_key[key2])*perm*perm_key[key1]).tuple())
		else:
			tet.detach(TwoSubsimplices[c1])
			new_tet2.attach(F0,voisin,(perm*perm_key[key1]).tuple())
	if other_tet.Neighbor[TwoSubsimplices[c2]] != None:
		perm = other_tet.Gluing[TwoSubsimplices[c2]]
		voisin = other_tet.Neighbor[TwoSubsimplices[c2]]
		key1 = (other_tet,c2)
		if voisin == tet or voisin == other_tet:
			key2 = (voisin,FaceIndex[perm.image(TwoSubsimplices[c2])])
			new_tet2.attach(F3,tet_key[key2],(inv(perm_key[key2])*perm*perm_key[key1]).tuple())
		else:
			other_tet.detach(TwoSubsimplices[c2])
			new_tet2.attach(F3,voisin,(perm*perm_key[key1]).tuple())
	new_tet0.fill_edge_params(z0*w0/(z0 - ComplexSquareRootCombination.One() + w0))
	new_tet1.fill_edge_params((z0 - ComplexSquareRootCombination.One() + w0)/w0)
	new_tet2.fill_edge_params((w0 - ComplexSquareRootCombination.One())/(z0 - ComplexSquareRootCombination.One() + w0))
	# Now we have to think about if the face was glued to itself, or if there was rotation.
	if face_glued_to_self and not face_rotation:
		# vert0 is the vertex fixed by the face gluing map
		if vert0 == b1:	
			new_tet0.Symmetries = [Perm4((0,1,2,3)),Perm4((3,2,1,0))]
			new_tet0.detach(F3)
			new_tet1.detach(F3)
			new_tet1.attach(F3,new_tet1,[2,1,0,3])
			# We're going to get rid of new_tet2, but some other faces of new_tet0 or new_tet1 could be
			# attached to new_tet2. In that case we detach, then glue to a corresponding face of
			# new_tet0 or new_tet1. 
			if new_tet0.Neighbor[F1] == new_tet2:
				perm = new_tet0.Gluing[F1]
				new_tet0.detach(F1)
				new_tet0.attach(F1,new_tet1,(Perm4((2,1,3,0))*perm).tuple())
			if new_tet1.Neighbor[F0] == new_tet2:
				perm = new_tet1.Gluing[F0]
				new_tet1.detach(F0)
				new_tet1.attach(F0,new_tet1,Perm4((2,1,3,0))*perm)

			new_tets = [new_tet0,new_tet1]
			for T in triang:
				if T != tet:
					new_tets.append(T)
			index_tets(new_tets)
			return new_tets
		elif vert0 == a1:
			new_tet1.Symmetries = [Perm4((0,1,2,3)),Perm4((2,3,0,1))]
			new_tet1.detach(F3)
			new_tet0.detach(F3)
			new_tet0.attach(F3,new_tet0,[0,2,1,3])
			# We're going to get rid of new_tet2, but some other faces of new_tet0 or new_tet1 could be
			# attached to new_tet2. In that case we detach, then glue to a corresponding face of
			# new_tet0 or new_tet1.
			if new_tet0.Neighbor[F1] == new_tet2:
				perm = new_tet0.Gluing[F1]
				new_tet0.detach(F1)
				new_tet0.attach(F1,new_tet0,(Perm4((2,3,0,1))*perm).tuple())
			if new_tet1.Neighbor[F0] == new_tet2:
				perm = new_tet1.Gluing[F0]
				new_tet1.detach[F0]
				new_tet1.attach(F0,new_tet0,(Perm4((2,3,0,1))*perm).tuple())

			new_tets = [new_tet0,new_tet1]
			for T in triang:
				if T != tet:
					new_tets.append(T)
			index_tets(new_tets)
			return new_tets
		elif vert0 == c1:
			new_tet2.Symmetries = [Perm4((0,1,2,3)),Perm4((3,2,1,0))]
			new_tet2.detach(F2)
			new_tet0.detach(F0)
			new_tet0.attach(F0,new_tet0,[0,2,1,3])
			# We're going to get rid of new_tet1, but some other faces of new_tet0 or new_tet2 could be
			# attached to new_tet1. In that case we detach, then glue to a corresponding face of
			# new_tet0 or new_tet2.
			if new_tet0.Neighbor[F1] == new_tet1:
				perm = new_tet0.Gluing[F1]
				new_tet0.detach(F1)
				new_tet0.attach(F1,new_tet0,(Perm4((2,0,1,3))*perm).tuple())
			if new_tet2.Neighbor[F0] == new_tet1:
				perm = new_tet2.Gluing[F0]
				new_tet2.detach(F0)
				new_tet2.attach(F0,new_tet0,(Perm4((2,0,1,3))*perm).tuple())

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
		# Now it could be that F1 and F2 are glued to new_tet1 or new_tet2,
		# in that case need to detach then attach to correct face of new_tet0.... to be done later.
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





