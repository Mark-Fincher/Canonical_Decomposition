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

"""
Need interval arithmetic for this.

def check_2_to_3_possible(self,tet,face):
	for edge in OneSubsimplices:
		if is_subset(edge,face):
			z = tet.edge_params[edge]
			w = tet.Neighbor[face].edge_params[tet.Gluing[face].image(edge)]
			if (z*w).imag ...

"""


def two_to_three(triang,tet,face):
	#currently triang should be a list of tets, not a CuspedOrbifold object. That could change.
	custom_PickAnEdge = {F0:(3,1),F1:(2,0),F2:(0,1),F3:(2,1)}
	a1,b1 = custom_PickAnEdge[face]
	d1 = FaceIndex[face]
	c1 = ({0,1,2,3}-{a1,b1,d1}).pop()
	z0 = tet.edge_params[bitmap((a1,b1))]
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
	new_tets = [new_tet0,new_tet1,new_tet2]
	for T in triang:
		if T != tet and T != other_tet:
			new_tets.append(T)
	index_tets(new_tets)
	return new_tets





