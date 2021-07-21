# File for testing and learning related to the canonical decomposition project.
"""
from tetrahedron import*
from perm4 import*
from Exact_Arithmetic import*
from simplex import*
from Dest_to_Triang import*
"""
from The_Algorithm import*

"""
z = ComplexSquareRootCombination(SquareRootCombination([(1,Fraction(1/2))]),SquareRootCombination([(Fraction(3,1),Fraction(1/2))]))
tet0 = Tetrahedron()
tet0.fill_edge_params(z)
tet0.Index = 0
#tet0.Symmetries = [Perm4((0,1,2,3)),Perm4((0,2,3,1)),Perm4((0,3,1,2))]
tet1 = Tetrahedron()
tet1.fill_edge_params(z)
tet1.Index = 1
tet1.Symmetries = tet0.Symmetries
tet0.attach(F0,tet0,[0,1,3,2])
show_triangulation([tet0])
New_triang = two_to_three([tet0],tet0,F0)
show_triangulation(New_triang)
"""



#Dest = [0,0,0,0]
#Dest = [0,1,1,0,1,0,0,2,3,2,2,1,2,3,3,3]
#Dest = [0,1,2,3,2,2,0,2,1,0,1,1,4,3,3,0,3,4,4,4]
Dest = [0,1,2,1,2,3,0,0,1,0,4,2,4,5,1,4,3,2,5,3,5,4,3,6,7,6,6,5,6,7,7,7]

tets_list = full_snappy_triang(Dest)

#show_triangulation(tets_list)

next_list = two_to_three(tets_list,tets_list[0],F1)

#show_triangulation(next_list)

tet0 = CuspedOrbifold(next_list).Tetrahedra[0]

#print(check_2_to_3_possible([tet0],tet0,F1))

next_list = two_to_three([tet0],tet0,F1)

show_triangulation(next_list)


"""
print(tet0.tilt(V0).evaluate() < 0)
print(tet0.tilt(V1).evaluate() < 0)
print(tet0.tilt(V2).evaluate() < 0)
print(tet0.tilt(V3).evaluate() < 0)
"""



#print(check_2_to_3_possible(tets_list,tets_list[0],F1))



"""
z0 = ComplexSquareRootCombination(SquareRootCombination([(1, Fraction(1, 2))]), SquareRootCombination([(Fraction(3,1), Fraction(1, 2))]))
w0 = z0
print(z0*w0/(z0 - ComplexSquareRootCombination.One() + w0))
print(z0*w0)
print(z0 - ComplexSquareRootCombination.One() + w0)
print(ComplexSquareRootCombination.One()/ComplexSquareRootCombination(SquareRootCombination([]), SquareRootCombination([(Fraction(3,1), Fraction(1, 1))])))
"""
