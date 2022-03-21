"""
Isomorpism signature for simplicial orbifolds.
"""

from tetrahedron import*
from arrow import*
from edge import*
from vertex import*
from simplex import*
from corner import*
from perm4 import*
from tetrahedron import*
from SimplicialOrbifold import*

"""
Change a CuspedOrbifold object into a SimplicialOrbifold object. On a high level, this just means
forget about the geometry. We take the isomorphism signature of SimplicialOrbifold objects,
not CuspedOrbifold objects, which is why we might want to make this change. 
"""
def hyperbolic_to_simplicial(cusped_orb):
	tets_list = cusped_orb.Tetrahedra
	for tet in tets_list:
		tet.remove_extra_glued_to_self()
		tet.edge_params = {E01:None,E23:None,E02:None,E13:None,E03:None,E12:None}
		tet.horotriangles = {V0:None, V1:None, V2:None, V3:None}
		for one_subsimplex in OneSubsimplices:
			tet.edge_labels[one_subsimplex] = tet.Class[one_subsimplex].LocusOrder
		tet.clear_Class()
	return SimplicialOrbifold(tets_list)

"""
We want to find all canonical relabelings of a SimplicialOrbifold object orb, then encode 
those relabelings as strings, then get the lexicographically smallest of those strings, which 
we'll call the isomorphism signature. This is Burton's isomorphism signature, but for simplicial
orbifolds instead of manifolds.

A SimplicialOrbifold object automatically has a labeling of tetrahedra and vertices, and the face
gluing maps are described in terms of the vertex labeling. If we change the labeling, it's clearly
the same simplicial orbifold. We can even change some face gluings using symmetries, and we still
have the same orbifold. To describe the kinds of labelings we're looking for, we must define
the destination sequence and the type sequence.

The DESTINATION SEQUENCE corresponding to a labeling of a simplicial orbifold (having n tetrahedra)
is the sequence

D = D_{0,0}, D_{0,1}, D_{0,2}, D_{0,3}, D_{1,0}, D_{1,1}, ... , D_{n-1,2}, D_{n-1,3}

where D_{t,f} is the index of the tetrahedron glued to face f of tetrahedron t, or it's None
if f is not glued to anything. The TYPE SEQUENCE of the labeling is the sequence

T = T_{0,0}, T_{0,1}, T_{0,2}, T_{0,3}, T_{1,0}, T_{1,1}, ... , T_{n-1,2}, T_{n-1,3}

where T_{t,f} = 0 if face f of tetrahedron t is glued to None, T_{t,f} = 1 if A_{t,f} = k
where k != 0 and A_{t,f} is the first instance of k in A, or T_{t,f} = 2 if it doesn't equal
0 or 1. 

We say a labeling is CANONICAL if the following are satisfied.

1. For 0 < i < j < n, the first instance of i in D appears before the first instance of j.
2. For each face f of each tetrahedron t, if T_{t,f} = 1 then face f of tetrahedron t is glued
to tetrahedron D_{t,f} via the identity permutation.
3. For each face f of each tetrahedron t, if T_{t,f} = 1 and tetrahedron t has a symmetry taking face f to 
a different face f', then f < f'.
4. For each face f of each tetrahedron t, if T_{t,f} = 2 and face f of tetrahedron t is not glued to a
face of type 1, and tetrahedron t has a symmetry taking face f to a different face f', then f < f'.

It can be seen that for each choice of a tetrahedon to be labeled 0 and choice of labeling of its vertices,
there is a corresponding canonical labeling, and every canonical labeling of a simplicial orbifold arises
in this way.

To keep track of relabelings, we do not create new SimplicialOrbifold objects. Instead, we create
two lists, tets and perms, where tets[i] is the element of self.Tetrahedra which is now labeld i, and
perms[i] represents  the relabeling of tets[i], i.e. the vertex j becomes perms[i][j]. Note that
in general tets[i].Index != i, i.e. the Index attribute is not updated.

We use three other lists to describe simpicial orbifold data, G gluing sequence, S symmetry sequence,
and E edge label sequence.

G = G_{0,0}, G_{0,1}, G_{0,2}, G_{0,3}, G_{1,0}, ... G_{n-1,2}, G_{n-1,3}

where G_{t,f} is the permutation gluing face f of tetrahedron t to whatever tetrahedron it's glued to.
We list all permutations in S4 lexicographically and take G_{t,f} to be the integer which is the index
of the corresponding permutation.

S = S_0, S_1, ... S_n-1

where S_t is the group of symmetries of tetrahedron t, encoded as integers based on an ordering
of groups given in the list ordered_subgrps_S4 below.

E = E_{0,0}, E_{0,1}, ... E_{0,5}, E_{1,0}, ... , E_{n-1,0}, E_{n-1,1}, ... E_{n-1,5}

where E_{t,e} is the edge label of edge e in tetrahedron t. We consider e an integer 0 <= e <= 5
corrseponding to the ordering E01, E02, E03, E12, E13, E23. This is listed in simplex.py.

In the following, it's often best to think of all implicit face gluings (due to symmetries) as actually
being explicit, until we fix the dest seq, type seq, and gluing seq. To that end we use the 
tetrahedron method true_glued. 	
"""

"""
For ordered lists of permutations, we will use _rawS4 and _rawA4 from perm.py.
The ordered list of subgroups of even permutations is as follows.
"""

_rawA4 = [ (0,1,2,3),
               (0,2,3,1),
               (0,3,1,2),
               (1,0,3,2),
               (1,2,0,3),
               (1,3,2,0),
               (2,0,1,3),
               (2,1,3,0),
               (2,3,0,1),
               (3,0,2,1),
               (3,1,0,2),
               (3,2,1,0)]

ordered_subgrps_S4 = (	[(0,1,2,3)],
						[(0,1,2,3),(1,0,3,2)],
						[(0,1,2,3),(2,3,0,1)],
						[(0,1,2,3),(3,2,1,0)],
						[(0,1,2,3),(0,2,3,1),(0,3,1,2)],
						[(0,1,2,3),(2,1,3,0),(3,1,0,3)],
						[(0,1,2,3),(1,3,2,0),(3,0,2,1)],
						[(0,1,2,3),(1,0,3,2),(2,3,0,1),(3,2,1,0)],
						_rawA4)

"""
Encode integers as printable characters. This is exactly as in the appendix to Burton's
"The Pachner graph and the simplification of 3-sphere triangulations".
"""
pi = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 
	'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
	'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 
	'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
	'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '-')

"""
def enc(i):
	# input: integer i.
	# output: string encoding i, using pi above.
	# First see what i is in base 64.
"""	


"""
Create labeled data for a canonical labeling.

Input:

orb, a SimplicialOrbifold object
tet0, a tetrahedron in orb.Tetrahedra which will be relabeled 0
perm0, a perm4 object representing the relabeling of vertices of tet0

Output:

The lists D,T,G as above, along with lits tets and perms.
"""
def get_tets_perms_D_T_G(orb,tet0,perm0):
	# Assuming tet0 should be given label 0 and perm0 should be the relabeling permutation
	# of tet0, find the corresponding canonical labeling in terms of the sequences
	# D,T,G,S,E. Lists which determine the relabeling are tets and perms. Meaning the tet
	# with label i after relabeling is tets[i], and vertex j of tets[i] becomes perms[i][j]
	# in the new labeling of vertices of tets[i]. We start only knowing the first elements
	# of tets and perms, and build the rest as we go along.
	tets = [tet0]
	perms = [perm0]
	D = 4*len(tets)*[None]
	T = 4*len(tets)*[None]
	G = 4*len(tets)*[None]
	for i in range(len(orb.Tetrahedra)):
		for j in range(4):
			# We are interested in face j of tet i, w.r.t. the new labeling.
			k = perms[i].inverse()[j]
			# In the original labeling, face j is face k.
			Neighbor, Gluing = tets[i].true_glued(TwoSubsimplices[k])
			if Neighbor in tets:
				if T[4*i + j] == 2:
					# Then, since the type was already assigned, this face must be glued to a
					# face of type 1. In which case we don't adjust any gluings even though
					# there might be symmetries. And note the face gluing map must be the identity.
					D[4*i + j] = tets.index(Neighbor)
					G[4*i + j] = 0
				else:
					# In this case the face is not type 1, nor is it glued to a type 1 face.
					for sym in tets[i].Symmetries:
						if (perms[i]*sym*perms[i].inverse())[j] < j:
							# Then face j is taken to some other face of lesser index by a symmetry.
							# So the gluing of face j should be implicit. 
							D[4*i + j] = None
							T[4*i + j] = 0
							G[4*i + j] = None
							break
					else:
						# There is no symmetry taking face j to a face of lesser index. It could be there
						# is no symmetry taking face j to a different face at all. In either case, we
						# have an explicit face gluing.
						m = tets.index(Neighbor)
						D[4*i + j] = m
						T[4*i + j] = 2
						G[4*i + j] = _rawS4.index((perms[m]*Gluing*perms[i].inverse()).tuple())
			elif Neighbor is not None and Neighbor not in tets:
				# Then face j of tet i is type 1.
				tets.append(Neighbor)
				# The new permutation, which describe the relabeling of Neighbor, must satisfy
				# new_perm*Gluing*perms[i].inverse() == identity for the labeling to be canonical.
				new_perm = perms[i]*Gluing.inverse()
				perms.append(new_perm)
				D[4*i + j] = tets.index(Neighbor)
				G[4*i + j] = 0 
				T[4*i + j] = 1
				# Set the type of the face f glued to face j of tet i to 2. This helps
				# us know later that the gluing of f should remain explicit, even if
				# there's a symmetry taking it to a face with lesser index. Note that,
				# since face j of tet i is glued with the identity map, the index of the
				# face it's glued to is also j. 
				T[4*tets.index(Neighbor) + j] = 2
			else:
				# In this case it really is a boundary face.
				D[4*i + j] = None
				T[4*i + j] = 0
				G[4*i + j] = None
	return (tets,perms,D,T,G)


def get_S_E(tets,perms):
	# Now let's get the set of symmetries, S, and edge labels, E. We only need one edge label
	# from each edge class.
	S = []
	E = []
	seen_edge_classes = []
	for i in range(len(tets)):
		sym_group = [(perms[i]*sym*perms[i].inverse()).tuple() for sym in tets[i].Symmetries]
		S.append(ordered_subgrps_S4.index(sorted(sym_group)))
		for edge in OrderedEdges:
			edge_class = tets[i].Class[perms[i].inverse().image(edge)]
			if edge_class not in seen_edge_classes:
				seen_edge_classes.append(edge_class)
				E.append(tets[i].edge_labels[perms[i].inverse().image(edge)])
	return (S,E)


"""
Assuming that we have D, T, G, we want to remove redundant information.
S and E are already as concise as possible, so we don't do anything with
them here.

Input: lists of integers D,T,G,S,E.

Output: D,T,G, with some elements removed.
"""
def shorten_lists(D,T,G,S):
	# The indices of the elements of D to be removed.
	faces_to_remove = []
	for i in len(S):
		for j in range(4):
			if D[4*i + j] is not None:
				# Then this face is glued to something. Let's check if it's glued to a
				# (tet,face) pair which is lexicographically smaller then (i,j).
				if 4*D[4*i + j] + _rawS4[G[4*i + j]][j] < 4*i + j:
					faces_to_remove.append(4*i + j)
					continue
			# Now we check if this face is identified by a symmetry with a smaller face.
			sym_group = ordered_subgrps_S4[S[i]]
			for sym in sym_group:
				if sym[j] < j:
					faces_to_remove.append(4*i + j)
					break
	new_D = [D[k] for k in len(D) if k not in faces_to_remove]
	new_T = [T[k] for k in len(D) if k not in faces_to_remove]
	new_G = [G[k] for k in len(D) if k not in faces_to_remove]
	return (new_D,new_T,new_G)
