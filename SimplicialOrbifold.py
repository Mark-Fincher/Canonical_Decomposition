from tetrahedron import*
from arrow import*
from edge import*
from vertex import*
from simplex import*
from corner import*
from perm4 import*


class SimplicialOrbifold:
	def __init__(self,tets_list):
		self.Tetrahedra = tets_list
		for i in range(len(self.Tetrahedra)):
			self.Tetrahedra[i].Index = i
		self.Edges = []
		self.Faces = []
		self.Vertices = []
		#self.build_vertex_classes()
		#self.build_edge_classes()


	"""
	We create edge classes just like we do for a CuspedOrbifold. The only difference is
	that we don't have to figure out the LocusOrder at the end using the geometry.
	"""
	def build_edge_classes(self):
		for tet in self.Tetrahedra:
			for one_subsimplex in OneSubsimplices:
				if tet.Class[one_subsimplex] is None:
					newEdge = Edge()
					self.Edges.append(newEdge)
					a = Arrow(one_subsimplex,RightFace[one_subsimplex],tet)
					first_arrow = a.copy()
					sanity_check = 0
					unfinished_walk = True
					while unfinished_walk:
						if sanity_check > 6*len(self.Tetrahedra):
							print('Bad gluing data: could not construct edge link.')
						newEdge.Corners.append(Corner(a.Tetrahedron, a.Edge))
						if a.true_next() is None:
							print('Hit boundary. Did not construct edge link.')
							break
						else:
							#Check if a is now first_arrow, or the image of a under 
							#a symmetry is first_arrow.
							if a.Tetrahedron is first_arrow.Tetrahedron:
								for sym in a.Tetrahedron.Symmetries:
									if ( sym.image(a.Edge) == first_arrow.Edge 
										and 
										sym.image(a.Face) == first_arrow.Face):
										#Then we've finished walking around the edge,
										#and we want to stop the while loop.
										unfinished_walk = False
										break
							sanity_check = sanity_check + 1
					#We've finished walking around the edge, and now we want to assign newEdge to all
					#edges in this class. This might not just be the edges in newEdge.Corners; have
					#to apply symmetries to get other edges in the class.
					for corner in newEdge.Corners:
						for sym in corner.Tetrahedron:
							corner.Tetrahedron.Class[sym.image(corner.Subsimplex)] = newEdge
		for i in range(len(self.Edges)):
			self.Edges[i].Index = i


	"""
	Isomorphism signature.

	Briefly: we want to find all canonical relabelings, then encode those relabelings as strings, then
	get the lexicographically smallest of those strings, which we'll call the isomorphism signature.

	A SimplicialOrbifold object automatically has a labeling of tetrahedra and vertices, and the face
	gluing maps are described in terms of the vertex labeling. If we change the labeling, it's clearly
	the same simplicial orbifold. To describe the kinds of labelings we're looking for, we must define
	the destination sequence and the type sequencce.

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

	1. For 0 < i < j < n, the first instance of i in A appears before the first instance of j.
	2. For each face f of each tetrahedron t, if T_{t,f} = 1 then face f of tetrahedron t is glued
	to tetrahedron A_{t,f} via the identity permutation.
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
	in general tets[i].Index != i.

	We use three other lists to describe simpicial orbifold data, G gluing sequence, S symmetry sequence,
	and E edge label sequence.

	G = G_{0,0}, G_{0,1}, G_{0,2}, G_{0,3}, G_{1,0}, ... G_{n-1,2}, G_{n-1,3}

	where G_{t,f} is the permutation gluing face f of tetrahedron t to whatever tetrahedron it's glued to.
	We list all permutations in S4 lexicographically and take G_{t,f} to be the integer which is the index
	of the corresponding permutation.

	S = S_0, S_1, ... S_n-1

	where S_t is the group of symmetries of tetrahedron t, encoded as integers like the elements of G.

	E = E_{0,0}, E_{0,1}, ... E_{0,5}, E_{1,0}, ... , E_{n-1,0}, E_{n-1,1}, ... E_{n-1,5}

	where E_{t,e} is the edge label of edge e in tetrahedron t. We consider e an integer 0 <= e <= 5
	corrseponding to the ordering E01, E02, E03, E12, E13, E23.

	In the below, it's often best to think of all implicit face gluings (due to symmetries) as actually
	being explicit, until we fix the dest seq, type seq, and gluing seq. To that end we use the 
	tetrahedron method true_glued. 	
	"""

	def get_relabeling(self,tet0,perm0):
		#Assuming tet0 should be given label 0 and perm0 should be the relabeling permutation
		#of tet0, find the corresponding canonical labeling in terms of the sequences
		#D,T,G,S,E. Lists which determine the relabeling are tets and perms. Meaning the tet
		#with label i after relabeling is tets[i], and vertex j of tets[i] becomes perms[i][j]
		#in the new labeling of vertices of tets[i].
		tets = [tet0]
		perms = [perm0]
		D = 4*len(tets)*[None]
		T = 4*len(tets)*[None]
		G = 4*len(tets)*[None]
		for i in range(len(self.Tetrahedra)):
			for j in range(4):
				#We are interested in face j of tet i, w.r.t. the new labeling.
				k = perms[i].inverse().image(j)
				#In the original labeling, face j is face k.
				Neighbor, Gluing = tets[i].true_glued(TwoSubsimplices[k])
				if Neighbor in tets:
					if T[4*i + j] == 2:
						#Then, since the type was already assigned, this face must be glued to a
						#face of type 1. In which case we don't adjust any gluings even though
						#there might be symmetries.
						D[4*i + j] = tets.index(Neighbor)
						G[4*i + j] = Gluing*perms[i].inverse()
					else:
						#In this case the face is not type 1, nor is it glued to a type 1 face.
						for sym in tets[i].Symmetries:
							if (perms[i]*sym*perms[i].inverse())[j] < j:
								#Then face j is taken to some other face of lesser index by a symmetry.
								#So the gluing of face j should be implicit. 
								D[4*i + j] = None
								T[4*i + j] = 0
								G[4*i + j] = None
								break
						else:
							#There is no symmetry taking face j to a face of lesser index. It could be there
							#is no symmetry taking face j to a different face at all. In either case, we
							#have an explicit face gluing.
							D[4*i + j] = tets.index(Neighbor)
							T[4*i + j] = 2
							G[4*i + j] = Gluing*perms[i].inverse()
				elif Neighbor is not None and Neighbor is not in tets:
					#Then face j of tet i is type 1.
					tets.append(Neighbor)
					perms.append(Gluing*perms[i].inverse())
					D[4*i + j] = tets.index(Neighbor)
					G[4*i + j] = Gluing*perms[i].inverse() 
					T[4*i + j] = 1
					#Set the type of the face f glued to face j of tet i to 2. This helps
					#us know later that the gluing of f should remain explicit, even if
					#there's a symmetry taking it to a face with lesser index. 
					T[4*tets.index(Neighbor) + j] = 2
				else:
					#In this case it really is a boundary face.
					D[4*i + j] = None
					T[4*i + j] = 0
					G[4*i + j] = None
		#Now let's get the set of symmetries, S, and edge labels, E.
		S = []
		E = []
		for i in range(len(tets)):
			S.append([])
			for sym in tets[i].Symmetries:
				S[i].append(perms[i]*sym*perms[i].inverse())
			for edge in OrderedEdges:
				E.append(tets[i].edge_group_labels[perms[i].inverse().image(edge)])
		return (D,T,G,S,E)


	def encode_as_string(self, D,T,G,S,E):
		#The lists D,T,G,S,E as they are returned from get_relabeling determine
		#a relabeled simplicial orbifold. We can encode this data more compactly
		#as a string in the following way. First we clean them up somewhat.
		#Remove redundant information. If face f of tet t is glued to face f' of tet t'
		#and (t,f) < (t',f') lexicographically, we don't need to include face gluing
		#data of (t',f') because we will already see it in the gluing data of (t,f).
		#So we remove D[t',f'] from the destination sequence, then we remove the
		#corresponding term from the type sequence and the gluing sequence.
		#There is also redundant information in the edge label sequence E. If edge e of tetrahedron
		#t is mapped to edge e' of tetrahedron t' via a symmetry or face gluing map, we only need
		#to label the lexicographically smallest of the two edges, since the labeling of the other
		#can be inferred. 
		
		#There is redundant information in the edge label sequence E. If edge e of tetrahedron
		#t is mapped to edge e' of tetrahedron t' via a symmetry or face gluing map, we only need
		#to label the lexicographically smallest of the two edges, since the labeling of the other
		#can be inferred. 
		to_remove = []
		for i in range(len(tets)):
			for j in range(6):
				edge = OrderedEdges[j]
				k = TwoSubsimplices.index(RightFace[edge])
				if D[4*i + k] is not None:
					if D[4*i + k] < i:
						to_remove.append(4*i + j)
					elif D[4*i + k] == i:
						


		to_remove = []
		for i in range(len(tets)):
			for j in range(4):
				k = D[4*i + j]
				if k is not None:
					l = G[4*i + j].image(j)
					if 4*k + l < 4*i + j:
						to_remove.append(4*i + j)
		for k in to_remove:
			D.remove(D[k])
			T.remove(T[k])
			G.remove(G[k])
		#Now we change G and S so that each permutation is represented by an integer, its index
		#in the lexicographically ordered list of all elements of S4. In perm4.py this list
		#is called _rawS4.
		G = [_rawS4.index(perm.tuple()) for perm in G]
		for i in len(S):
			int_S_i = [_rawS4.index(perm_tuple()) for perm in S[i]]
			S[i] = sorted(int_S_i)


