"""
This program takes a destination sequence, turns it into a "naive ideal triangulation", then turns that into a collection 
of objects in the snappy tetrahedron class, which more verbosely describes the triangulation.
"""

from perm4 import *
from tetrahedron import *
from simplex import *
from Exact_Arithmetic import *
from CuspedOrbifold import*

def fill_face(f,i):
	# If the face f has a block filled in at f[i], then this function will try fill in the other entries with blocks
	# until either the whole face is filled, or you reach another entry which already had a block.

	if i != 5 and f[i+1] == 0:
		if f[i][1] == 0:
			f[i+1] = [D[4*f[i][0]+3],1]
		if f[i][1] == 1:
			f[i+1] = [D[4*f[i][0]],0]
		fill_face(f,i+1)
	if i == 5 and f[0] == 0:
		if f[5][1] == 0:
			f[0] = [D[4*f[5][0]+3],1]
		if f[5][1] == 1:
			f[0] = [D[4*f[5][0]],0]
		fill_face(f,0)

def filled_tet(face_label,face):
	# given that face_label of a tet is filled in partially according to face,
	# return the completely filled in tet this determines
	filled = [0,0,0,0]
	for i in range(6):
		if face[i] != 0:
			fill_face(face,i)
	# Now that particular face of the tet is filled, and we fill the others based on that
	# and what its label is
	if face_label == 0:
		filled[0] = face
		filled[1] = [0,[D[4*face[0][0]+1],1],0,0,0,0]
		filled[2] = [[D[4*face[1][0]+2],0],0,0,0,0,0]
		filled[3] = [0,[D[4*face[4][0]+1],1],0,0,0,0]
		fill_face(filled[1],1)
		fill_face(filled[2],0)
		fill_face(filled[3],1)
	if face_label == 1:
		filled[1] = face
		filled[2] = [0,[D[4*face[0][0]+1],1],0,0,0,0]
		filled[0] = [[D[4*face[1][0]+2],0],0,0,0,0,0]
		filled[3] = [[D[4*face[3][0]+2],0],0,0,0,0,0]
		fill_face(filled[0],0)
		fill_face(filled[2],1)
		fill_face(filled[3],0)
	if face_label == 2:
		filled[2] = face
		filled[0] = [0,[D[4*face[0][0]+1],1],0,0,0,0]
		filled[1] = [[D[4*face[1][0]+2],0],0,0,0,0,0]
		fill_face(filled[0],1)
		fill_face(filled[1],0)
		filled[3] = [0,[D[4*filled[0][4][0]+1],1],0,0,0,0]
		fill_face(filled[3],1)
	if face_label == 3:
		filled[3] = face
		filled[0] = [0,0,0,0,[D[4*face[1][0]+2],0],0]
		fill_face(filled[0],4)
		filled[1] = [0,[D[4*filled[0][0][0]+1],1],0,0,0,0]
		filled[2] = [[D[4*filled[0][1][0]+2],0],0,0,0,0,0]
		fill_face(filled[1],1)
		fill_face(filled[2],0)
	return filled

def update_blocks(triang):
	for face in triang[-1]:
		for block in face:
			if block not in seen_blocks:
				seen_blocks.append(block)

def add_tets(tet,triang):
	# Add tetrahedra around tet, unless it would be a copy of another tetrahedron already in triang.
	for i in range(4):
		if tet[i][0][1] == 0 and [tet[i][0][0],1] not in seen_blocks:
			triang.append(filled_tet(0,[0,[tet[i][0][0],1],0,0,0,0]))
			update_blocks(triang)
			add_tets(triang[-1],triang)
		if tet[i][0][1] == 1 and [tet[i][0][0],0] not in seen_blocks:
			triang.append(filled_tet(0,[0,[tet[i][0][0],0],0,0,0,0]))
			update_blocks(triang)
			add_tets(triang[-1],triang)



def dest_to_naive_triang(Dest):
	global D
	D = Dest
	T = [filled_tet(0,[[0,0],0,0,0,0,0])]
	global seen_blocks
	seen_blocks = []
	update_blocks(T)
	prev_length = 0
	while len(T) != prev_length:
		prev_length = len(T)
		for tet in T:
			add_tets(tet,T)
	return T

"""
Now, given a naive ideal triangulation, triang, we need to start turning it into a snappy triangulation. The first step
is the following function, find_spouse(block_position,triang). 

Recall that a block is a list of the form [i,0] or [i,1],
where i ranges through some allotted amount. We say that [i,0] and [i,1] are "spouses". A face f0
will be glued to a face f1 if and only if one contains the spouses of the other, and the locations of the spouses in the face
determine the particular face gluing map. 

A "block position" of a block is a triple (i,j,k) such that block == triang[i][j][k]. Since a block can occur in multiple
places in a tetrahedron, it can have multiple block positions.

find_spouse(block_position, triang) finds the block position of the spouse of the block at block_position in triang. If
the spouse occurs at more than one place, then it just returns the first place it finds it. We check in the same face first.  
"""

def find_spouse(block_position,triang):
	t,f,b = block_position
	block = triang[t][f][b]
	if block[1] == 0:
		spouse = [block[0],1]
	if block[1] == 1:
		spouse = [block[0],0]
	#First check in the same face as block
	for i in range(6):
		if triang[t][f][i] == spouse:
			return (t,f,i)
	#Now search in all the tets
	for i in range(len(triang)):
		for j in range(4):
			for k in range(6):
				if triang[i][j][k] == spouse:
					return (i,j,k)

"""
Now we define a function find_perm(block_position,spouse_position). For a block at block_position and a spouse at spouse_position,
we want to glue together the two faces they lie in, with the correct gluing map. This function returns the gluing map, as a permutation.
This is with respect to the labeling I was using on my tetrahedra, which is different from the labeling of snappy tetrahedra. Will have
to account for this later.
"""

def find_perm(block_position,spouse_position):
	a,b,c = block_position
	l,m,n = spouse_position
	# This function assumes c is actually 0! Don't input a block_position with c != 0.
	# We will only use b, m, and n. Using the following permutations allows you to be more compact.
	e = Perm4((0,1,2,3))
	s1 = Perm4((1,2,0,3))
	s2 = Perm4((2,0,1,3))
	s3 = Perm4((3,1,0,2))
	t1 = Perm4((0,1,3,2))
	t2 = Perm4((0,3,2,1))
	t3 = Perm4((0,2,1,3))
	s = [e,s1,s2,s3]
	if n == 1:
		perm = s[m]*t3*inv(s[b])
		return perm 
	elif n == 3:
		perm = s[m]*t2*inv(s[b])
		return perm
	elif n == 5:
		perm = s[m]*t1*inv(s[b])
		return perm
	else:
		print('error in find_perm')
		return None

"""
Now it's time to start creating some snappy tetrahedra.
"""
def snappy_tets(triang):
	# sigma used to translate between the labeling of a niave tet and a snappy tet
	sigma = Perm4((2,3,0,1))
	tets = [Tetrahedron() for i in range(len(triang))]
	for i in range(len(tets)):
		tets[i].Index = i
	for i in range(len(triang)):
		for j in range(4):
			if tets[i].Neighbor[TwoSubsimplices[j]] == None:
				spouse_pos = find_spouse((i,sigma[j],0), triang)
				if tets[spouse_pos[0]].Neighbor[TwoSubsimplices[sigma[spouse_pos[1]]]] == None:
					permy = sigma*find_perm((i,sigma[j],0),spouse_pos)*sigma
					tets[i].attach(TwoSubsimplices[j],tets[spouse_pos[0]],(permy[0],permy[1],permy[2],permy[3]))
	return tets

def index_tets(tets_list):
	# add indices to tets
	index_count = 0
	for tet in tets_list:
		tet.Index = index_count
		index_count += 1


"""
Now we want to get the symmetries of each tet. At that point we'll be done encoding the triangulation. 
"""


def full_snappy_triang(Dest):
	triang = dest_to_naive_triang(Dest)
	tets = snappy_tets(triang)
	sigma = Perm4((2,3,0,1))
	e = Perm4((0,1,2,3))
	s1 = Perm4((1,2,0,3))
	s2 = Perm4((2,0,1,3))
	s3 = Perm4((3,1,0,2))
	t1 = Perm4((0,2,3,1))
	t2 = Perm4((0,3,1,2))
	s = [e,s1,s2,s3]
	for i in range(len(tets)):
		perms = []
		block_0 = triang[i][0][0]
		for j in range(4):
			for k in range(6):
				if triang[i][j][k] == block_0:
					if k == 0:
						perms.append(sigma*s[j]*sigma)
					if k == 2:
						perms.append(sigma*s[j]*t1*sigma)
					if k == 4:
						perms.append(sigma*s[j]*t2*sigma)
		tets[i].Symmetries = perms
		if len(perms) == 0:
			tets[i].Symmetries = [e]
	# Now add the shape parameters. Since the tetrahedra are regular, all edges get the same param, 1/2 + sqrt(3)*i/2
	for tet in tets:
		for i in range(6):
			tet.edge_params[OneSubsimplices[i]] = ComplexSquareRootCombination(SquareRootCombination([(1,Fraction(1/2))]),SquareRootCombination([(Fraction(3,1),Fraction(1/2))])) 
	# see the comments preceeding the function "unglue_some_faces" below.
	unglue_some_faces(tets)
	# add indices
	index_tets(tets)
	return tets

def dest_to_orb(Dest):
	# input a Dest Seq, returns the cusped orbifold object you get using Dest_to_Triang.py
	# and CuspedOrbifold.py.
	tets_list = full_snappy_triang(Dest)
	orb = CuspedOrbifold(tets_list)
	orb.DestSeq = Dest
	return orb


"""
If a tetrahedron has symmetries, it's not necessary to specify all the face gluings. This is because
you can use a symmetry to move an unglued face to a face which is glued, use that gluing data to say what
the unglued face should be glued to. In fact, we should exploit this. We should, in most cases, prefer to
have a face glued to None if possible, using symmetries. This makes things easier when we do 2-3 moves for
instance. The only exception should be if a face is glued to itself (this just doesn't cause problems).
The following function takes a list of tetrahedra and tries to unglue faces where the symmetries allow it.
"""

def unglue_some_faces(tets_list):
	for tet in tets_list:
		for face1 in TwoSubsimplices:
			for perm in tet.Symmetries:
				face2 = perm.image(face1)
				if face1 != face2 and (tet.Neighbor[face1] != tet or tet.Gluing[face1].image(face1) != face1):
					if tet.Neighbor[face1] != None and tet.Neighbor[face2] != None:
						tet.detach(face2)










