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
eye = ComplexSquareRootCombination(SquareRootCombination.Zero(),SquareRootCombination.One())
T1 = Tetrahedron()
T1.fill_edge_params(eye)
T2 = Tetrahedron()
T2.fill_edge_params(eye)
T1.attach(F0,T2,[0,1,3,2])
show_triangulation([T1,T2])
#New_triang = two_to_three(CuspedOrbifold([T1,T2]),T1,F0)
#show_triangulation(New_triang.Tetrahedra)
show_triangulation(two_to_three([T1,T2],T1,F0))
"""



#Dest = [0,0,0,0]
#Dest = [0,1,1,0,1,0,0,2,3,2,2,1,2,3,3,3]
#Dest = [0,1,2,3,2,2,0,2,1,0,1,1,4,3,3,0,3,4,4,4]
Dest = [0,1,2,1,2,3,0,0,1,0,4,2,4,5,1,4,3,2,5,3,5,4,3,6,7,6,6,5,6,7,7,7]

tets_list = full_snappy_triang(Dest)

show_triangulation(tets_list)

new_tets_list = two_to_three(tets_list,tets_list[0],bitmap((0,2,3)))

show_triangulation(new_tets_list)


"""
z0 = ComplexSquareRootCombination(SquareRootCombination([(1, Fraction(1, 2))]), SquareRootCombination([(Fraction(3,1), Fraction(1, 2))]))
w0 = z0
print(z0*w0/(z0 - ComplexSquareRootCombination.One() + w0))
print(z0*w0)
print(z0 - ComplexSquareRootCombination.One() + w0)
print(ComplexSquareRootCombination.One()/ComplexSquareRootCombination(SquareRootCombination([]), SquareRootCombination([(Fraction(3,1), Fraction(1, 1))])))
"""
