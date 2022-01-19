"""
This file contains the canonize function and a few other things which will be moved. The canonize
function will be updated soon.
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
proto_canonize tries to find a proto-canonical triangulation of an orbifold, which is analogous
to a proto-canonical triangulation of a manifold. If the orbifold's canonical decomposition is
a triangulation, as is usually the case, then its proto-canonical triangulation will actually
be its canonical decomposition

This function is written just like the proto_canonize function for manifolds in the SnapPy kernel, in 
canonize_part_1.c, with changes accounting for orbifolds.

This function will modify orb, so make a copy of it if you want to save it. If it succesfully makes
orb proto-canonical, it returns 1. Otherwise, it returns 0.
"""
def proto_canonize(orb):
	MAX_MOVES = 1000
	for i in range(MAX_MOVES):
		if attempt_cancellation(orb):
			print("cancelled flat tet(s)")
			continue
		if attempt_two_to_three(orb):
			print("did 2-3 move")
			continue
		if attempt_three_to_two(orb):
			print("did 3-2 move")
			continue
		if attempt_three_to_six(orb):
			print("did 3-6 move")
			continue
		if attempt_four_to_four(orb):
			print("did 4-4 move")
			continue
		if attempt_six_to_three(orb):
			print("did 6-3 move")
			continue
		if attempt_retriangulate_cube(orb):
			print("retriangulated a cube")
			continue
		#If none of the attempts work, then either it's proto-canonical, or the algorithm is stuck.
		#In either case, we break out of the loop.
		break
	if orb.is_proto_canonical():
		return 1
	else:
		return 0

def attempt_cancellation(orb):
	for edge in orb.Edges:
		#orb.cancel_tetrahedra(edge) will return 0 if "edge" doesn't belong to flat tet(s),
		#or if there are other reasons tet(s) around "edge" can't be cancelled. If they can
		#be cancelled, it will do so then return 1.
		if orb.cancel_tetrahedra(edge):
			return 1
	return 0

def attempt_two_to_three(orb):
	for tet in orb.Tetrahedra:
		for face in TwoSubsimplices:
			if concave_face(face,tet) and orb.check_two_to_three(face,tet):
				if orb.two_to_three(face,tet):
					return 1
				else:
					#This shouldn't happen
					print('error in attempt_two_to_three')
					return 0
	return 0

def attempt_three_to_two(orb):
	for edge in orb.Edges:
		if edge.valence() in {1,3} and concave_edge(edge):
			if orb.three_to_two(edge):
				return 1
	return 0

def attempt_three_to_six(orb):
	for tet in orb.Tetrahedra:
		for face in TwoSubsimplices:
			if concave_face(face,tet):
				if orb.three_to_six(face,tet):
					return 1
	return 0

def attempt_six_to_three(orb):
	for edge in orb.Edges:
		if concave_edge(edge) and orb.six_to_three(edge):
			return 1
	return 0

def attempt_four_to_four(orb):
	for edge in orb.Edges:
		if concave_edge(edge) and orb.four_to_four(edge):
			return 1
	return 0

def attempt_retriangulate_cube(orb):
	for tet in orb.Tetrahedra:
		for face in TwoSubsimplices:
			if concave_face(face,tet) and orb.retriangulate_cube(face,tet):
				return 1
	return 0


def concave_edge(edge):
	#Checks if there is any face adjacent to edge with positive tilt sum.
	#It should be that if any face adjacent to edge has positive tilt sum, then
	#they all do (except those which are glued to None), but it could be that,
	#for a given corner, both adjacent faces are glued to None. So we need to
	#loop through all the corners.
	#Will return 0 if any any tetrahedron adjacent to edge is flat, just because
	#we don't want to do a 3-2 move in that case.
	for corner in edge.Corners:
		tet = corner.Tetrahedron
		#First we look for flat tets. If any tetrahedron adjacent to edge is flat, 
		#we can't compute its tilts, so return 0.
		if tet.is_flat():
			return 0
	for corner in edge.Corners:
		tet = corner.Tetrahedron
		face = RightFace[corner.Subsimplex]
		for i in range(2):
			if tet.Neighbor[face] != None:
				other_tet = tet.Neighbor[face]
				other_face = tet.Gluing[face].image(face)
				if (tet.tilt(comp(face)) + other_tet.tilt(comp(other_face))).evaluate() > 0:
					return 1
			face = LeftFace[corner.Subsimplex]
	return 0

def concave_face(face,tet):
	#Check if face of tet is concave. Will return 0 if tet or tet.Neighbor[face]
	#is flat, just because we don't want to do a 2-3 move then.
	if tet.Neighbor[face] is None:
		return 0
	else:
		other_tet = tet.Neighbor[face]
		#If a tetrahedron is flat then we can't compute tilts, so return 0.
		if tet.is_flat() or other_tet.is_flat():
			return 0
		other_face = tet.Gluing[face].image(face)
		if (tet.tilt(comp(face)) + other_tet.tilt(comp(other_face))).evaluate() > 0:
			return 1
		else:
			return 0



# function that tries to find the canonical decomposition of the input_orb, for now only using 2-3 moves.
def old_canonize(input_orb):
	if input_orb.is_canonical is True:
		print('inputted orb is already canonical')
		return input_orb
	loop_limit = 100
	loop_count = 0
	active = [input_orb]
	while active:
		loop_count = loop_count + 1
		if loop_count < loop_limit:
			orb = active.pop()
			faces_seen = []
			for tet1 in orb.Tetrahedra:
				for face1 in TwoSubsimplices:
					if (tet1,face1) not in faces_seen and tet1.Neighbor[face1] != None:
						tet2,face2 = glued_to(tet1,face1)
						faces_seen.append((tet1,face1))
						faces_seen.append((tet2,face2))
						if (tet1.tilt(comp(face1)) + tet2.tilt(comp(face2))).evaluate() > 0:
							if check_2_to_3_possible(orb.Tetrahedra,tet1,face1) is True:
								#print("did 2-3")
								orb_copy = copy.deepcopy(orb)
								tet1_index = orb.Tetrahedra.index(tet1)
								new_orb = CuspedOrbifold(two_to_three(orb_copy.Tetrahedra,orb_copy.Tetrahedra[tet1_index],face1))
								new_orb.DestSeq = orb.DestSeq
								new_orb.PachnerPath = copy.deepcopy(orb.PachnerPath)
								new_orb.PachnerPath.append((orb,tet1,TwoSubsimplices.index(face1)))
								active = [new_orb] + active
								if new_orb.is_canonical is True:
									#print('Success! Found the canonical decomposition')
									return new_orb
		else:
			print('loop_limit reached')
			return
	print('failure: got stuck, could not reach the canonical decomposition with 2-3 moves')
	return


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
				new_tet1.attach(F0,new_tet1,(Perm4((2,1,3,0))*perm).tuple())
			
			# It could also be that some other tetrahedron in the triangulation is glued to
			# F0 of new_tet2. In that case we need to instead glue that tetrahedron to
			# F2 of new_tet1.
			voisin = new_tet2.Neighbor[F0]
			if voisin != None and voisin != new_tet0 and voisin != new_tet1:
				perm = new_tet2.Gluing[F0]
				new_tet2.detach(F0)
				if voisin != new_tet2:
					new_tet1.attach(F2,voisin,(perm*Perm4((3,1,0,2))).tuple())
				else:
					new_tet1.attach(F2,new_tet1,(inv(Perm4((3,1,0,2)))*perm*Perm4((3,1,0,2))).tuple())
			
			new_tets = [new_tet0,new_tet1]
			for T in triang:
				if T != tet:
					# clear classes and horotriangles! if they're not set to None then there are problems when you make
					# a cusped orbifold from this list of tets.
					T.clear_Class()
					T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
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

			# It could also be that some other tetrahedron in the triangulation is glued to
			# F0 of new_tet2. In that case we need to instead glue that tetrahedron to
			# F2 of new_tet0.
			voisin = new_tet2.Neighbor[F0]
			if voisin != None and voisin != new_tet0 and voisin != new_tet1:
				perm = new_tet2.Gluing[F0]
				new_tet2.detach(F0)
				if voisin != new_tet2:
					new_tet0.attach(F2,voisin,(perm*Perm4((2,3,0,1))).tuple())
				else:
					new_tet0.attach(F2,new_tet0,(inv(Perm4((2,3,0,1)))*perm*Perm4((2,3,0,1))).tuple())
			
			new_tets = [new_tet0,new_tet1]
			for T in triang:
				if T != tet:
					# clear classes and horotriangles! if they're not set to None then there are problems when you make
					# a cusped orbifold from this list of tets.
					T.clear_Class()
					T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
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

			# It could also be that some other tetrahedron in the triangulation is glued to
			# F0 of new_tet1. In that case we need to instead glue that tetrahedron to
			# F2 of new_tet0.
			voisin = new_tet1.Neighbor[F0]
			if voisin != None and voisin != new_tet0 and voisin != new_tet2:
				perm = new_tet1.Gluing[F0]
				new_tet1.detach(F0)
				if voisin != new_tet1:
					new_tet0.attach(F2,voisin,(perm*Perm4((1,2,0,3))).tuple())
				else:
					new_tet0.attach(F2,new_tet0,(inv(Perm4((1,2,0,3)))*perm*Perm4((1,2,0,3))).tuple())

			new_tets = [new_tet0,new_tet2]
			for T in triang:
				if T != tet:
					# clear classes and horotriangles! if they're not set to None then there are problems when you make
					# a cusped orbifold from this list of tets.
					T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
					T.clear_Class()
					new_tets.append(T)
			index_tets(new_tets)
			return new_tets
	elif face_rotation and not face_glued_to_self:
		new_tet0.detach(F3)
		new_tet0.detach(F0)
		new_tet0.attach(F0,new_tet0,[3,1,2,0])
		# Now it could be that F1 and F2 are glued to new_tet1 or new_tet2,
		# in that case need to detach then attach to appropriate face of new_tet0. Perm4((1,3,2,0))
		# is the rotation taking new_tet1 to new_tet0. Perm4((1,0,3,2)) is the rotation taking
		# new_tet2 to new_tet0.
		if new_tet0.Neighbor[F1] == new_tet1:
			perm = new_tet0.Gluing[F1]
			new_tet0.detach(F1)
			new_tet0.attach(F1,new_tet0,(Perm4((1,3,2,0))*perm).tuple())
		elif new_tet0.Neighbor[F1] == new_tet2:
			perm = new_tet0.Gluing[F1]
			new_tet0.detach(F1)
			new_tet0.attach(F1,new_tet0,(Perm4((1,0,3,2))*perm).tuple())
		if new_tet0.Neighbor[F2] == new_tet1:
			perm = new_tet0.Gluing[F2]
			new_tet0.detach(F2)
			new_tet0.attach(F2,new_tet0,(Perm4((1,3,2,0))*perm).tuple())
		elif new_tet0.Neighbor[F2] == new_tet2:
			perm = new_tet0.Gluing[F2]
			new_tet0.detach(F2)
			new_tet0.attach(F2,new_tet0,(Perm4((1,0,3,2))*perm).tuple())
		
		# Now consider the case that F1 and/or F2 are attached to None. This could happen,
		# since the rotational symmetries of the original two tetrahedra mean not all their faces
		# must be explicitly glued to something. And when we did the 2-3 move up to this point, we
		# didn't account for that.
		if new_tet0.Neighbor[F1] == None:
			# by rotating, this face maps to F0 of new_tet1, or, by rotating the other direction,
			# F0 of new_tet2. We see what those faces are glued to and attach F0 of new_tet0 to that.
			if new_tet1.Neighbor[F0] != None:
				voisin = new_tet1.Neighbor[F0]
				perm = new_tet1.Gluing[F0]
				new_tet1.detach(F0)
				if voisin == new_tet1:
					new_tet0.attach(F1,new_tet0,(inv(Perm4((3,0,2,1)))*perm*Perm4((3,0,2,1))).tuple())
				elif voisin == new_tet2:
					new_tet0.attach(F1,new_tet0,(inv(Perm4((1,0,3,2)))*perm*Perm4((3,0,2,1))).tuple())
				else:
					new_tet0.attach(F1,voisin,(perm*Perm4((3,0,2,1))).tuple())
			elif new_tet2.Neighbor[F0] != None:
				voisin = new_tet2.Neighbor[F0]
				perm = new_tet2.Gluing[F0]
				new_tet2.detach(F0)
				if voisin == new_tet1:
					new_tet0.attach(F1,new_tet0,(inv(Perm4((3,0,2,1)))*perm*Perm4((1,0,3,2))).tuple())
				elif voisin == new_tet2:
					new_tet0.attach(F1,new_tet0,(inv(Perm4((1,0,3,2)))*perm*Perm4((1,0,3,2))).tuple())
				else:
					new_tet0.attach(F1,voisin,(perm*Perm4((1,0,3,2))).tuple())
		if new_tet0.Neighbor[F2] == None:
			# by rotating, this face maps to F2 of new_tet1, or, in the other direction, F3 of new_tet2.
			# we now do the same thing as above.
			if new_tet1.Neighbor[F2] != None:
				voisin = new_tet1.Neighbor[F2]
				perm = new_tet1.Gluing[F2]
				new_tet1.detach(F2)
				if voisin == new_tet1:
					new_tet0.attach(F2,new_tet0,(inv(Perm4((3,0,2,1)))*perm*Perm4((3,0,2,1))).tuple())
				elif voisin == new_tet2:
					new_tet0.attach(F2,new_tet0,(inv(Perm4((1,0,3,2)))*perm*Perm4((3,0,2,1))).tuple())
				else:
					new_tet0.attach(F2,voisin,(perm*Perm4((3,0,2,1))).tuple())
			elif new_tet2.Neighbor[F3] != None:
				voisin = new_tet2.Neighbor[F3]
				perm = new_tet2.Gluing[F3]
				new_tet2.detach(F3)
				if voisin == new_tet1:
					new_tet0.attach(F2,new_tet0,(inv(Perm4((3,0,2,1)))*perm*Perm4((1,0,3,2))).tuple())
				elif voisin == new_tet2:
					new_tet0.attach(F2,new_tet0,(inv(Perm4((1,0,3,2)))*perm*Perm4((1,0,3,2))).tuple())
				else:
					new_tet0.attach(F2,voisin,(perm*Perm4((1,0,3,2))).tuple())
		# That concludes the case that F1 and/or F2 are attached to None.
		new_tets = [new_tet0]
		for T in triang:
			if T != tet and T != other_tet:
				# clear classes and horotriangles! if they're not set to None then there are problems when you make
				# a cusped orbifold from this list of tets.
				T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
				T.clear_Class()
				new_tets.append(T)
		index_tets(new_tets)
		return new_tets
	elif face_glued_to_self and face_rotation:
		new_tet0.detach(F3)
		new_tet0.detach(F0)
		new_tet0.attach(F0,new_tet0,[3,1,2,0])
		new_tet0.Symmetries = [Perm4((0,1,2,3)),Perm4((3,2,1,0))]
		# Again we have to consider the case that F1 and F2 of new_tet0 are glued to new_tet1
		# or new_tet2, since we're going to remove those two tets. Again we use the permutation
		# Perm4((1,3,2,0)), which is the rotation taking new_tet1 to new_tet0, and Perm4((1,0,3,2)),
		# which is the rotation taking new_tet2 to new_tet0.
		if new_tet0.Neighbor[F1] == new_tet1:
			perm = new_tet0.Gluing[F1]
			new_tet0.detach(F1)
			new_tet0.attach(F1,new_tet0,(Perm4((1,3,2,0))*perm).tuple())
		elif new_tet0.Neighbor[F1] == new_tet2:
			perm = new_tet0.Gluing[F1]
			new_tet0.detach(F1)
			new_tet0.attach(F1,new_tet0,(Perm4((1,0,3,2))*perm).tuple())
		# We don't worry about what F2 is attached to, since it's determined by F1 via the symmetry.
		# It could be that F1 is attached to None. In that case we want to rotate it to F0 of new_tet1 or
		# F0 of new_tet2 and see what they're attached to:
		if new_tet0.Neighbor[F1] == None:
			# by rotating, this face maps to F0 of new_tet1, or, by rotating the other direction,
			# F0 of new_tet2. We see what those faces are glued to and attach F0 of new_tet0 to that.
			if new_tet1.Neighbor[F0] != None:
				voisin = new_tet1.Neighbor[F0]
				perm = new_tet1.Gluing[F0]
				new_tet1.detach(F0)
				if voisin == new_tet1:
					new_tet0.attach(F1,new_tet0,(inv(Perm4((3,0,2,1)))*perm*Perm4((3,0,2,1))).tuple())
				elif voisin == new_tet2:
					new_tet0.attach(F1,new_tet0,(inv(Perm4((1,0,3,2)))*perm*Perm4((3,0,2,1))).tuple())
				else:
					new_tet0.attach(F1,voisin,(perm*Perm4((3,0,2,1))).tuple())
			elif new_tet2.Neighbor[F0] != None:
				voisin = new_tet2.Neighbor[F0]
				perm = new_tet2.Gluing[F0]
				new_tet2.detach(F0)
				if voisin == new_tet1:
					new_tet0.attach(F1,new_tet0,(inv(Perm4((3,0,2,1)))*perm*Perm4((1,0,3,2))).tuple())
				elif voisin == new_tet2:
					new_tet0.attach(F1,new_tet0,(inv(Perm4((1,0,3,2)))*perm*Perm4((1,0,3,2))).tuple())
				else:
					new_tet0.attach(F1,voisin,(perm*Perm4((1,0,3,2))).tuple())
			# That concludes the case that F0 of new_tet0 is attached to None.
		new_tets = [new_tet0]
		for T in triang:
			if T != tet:
				# clear classes and horotriangles! if they're not set to None then there are problems when you make
				# a cusped orbifold from this list of tets.
				T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
				T.clear_Class()
				new_tets.append(T)
		index_tets(new_tets)
		return new_tets
	else:
		new_tets = [new_tet0,new_tet1,new_tet2]
		for T in triang:
			if T != tet and T != other_tet:
				# clear classes and horotriangles! if they're not set to None then there are problems when you make
				# a cusped orbifold from this list of tets.
				T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
				T.clear_Class()
				new_tets.append(T)
		index_tets(new_tets)
		return new_tets





