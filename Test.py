# File for testing and learning related to the canonical decomposition project.
"""
from tetrahedron import*
from perm4 import*
from Exact_Arithmetic import*
from simplex import*
from Dest_to_Triang import*
"""
from The_Algorithm import*
import json


"""
Dest = [0,0,0,0]

tets_list = full_snappy_triang(Dest)

#show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

show_triangulation(orb.Tetrahedra)

tet0 = orb.Tetrahedra[0]

print(orb.isometries())
"""



"""
This is the 4 tet dest seq from Fact 5.5.
"""
"""
Dest = [0,1,1,0,1,0,0,2,3,2,2,1,2,3,3,3]

tets_list = full_snappy_triang(Dest)

orb = CuspedOrbifold(tets_list)

print("Dest seq (index 4) is")
print(Dest)
print(' ')
print("Regular triangulation is")
print(show_triangulation(orb.Tetrahedra))
print(' ')
print("Isometry group of regular triangulation is")
print(orb.isometries())
print(' ')
tet0 = orb.Tetrahedra[0]

print(tet0.tilt(V0))
print(tet0.tilt(V1))
print(tet0.tilt(V2))
print(tet0.tilt(V3))

# tet0.tilt(V2) is positive, and F2 is glued to itself, so that face is bad. Do 2-3 move through it.

#print(check_2_to_3_possible(orb.Tetrahedra,tet0,F2))

next_list = two_to_three(orb.Tetrahedra,orb.Tetrahedra[0],F2)

orb = CuspedOrbifold(next_list)
print("Canonical triangulation is")
print(show_triangulation(orb.Tetrahedra))

tet0 = orb.Tetrahedra[0]


print(tet0.tilt(V0))
print(tet0.tilt(V1))
print(tet0.tilt(V2))
print(tet0.tilt(V3))

# from checking tilt sums, this is canonical.
print(' ')
print("Isometry group of Canonical triangulation is")
print(orb.isometries())
print(' ')
"""


"""
This is one of the 8 tet dest seqs from Fact 5.5.
"""
"""
Dest = [0,1,2,1,2,3,0,0,1,0,4,2,4,5,1,4,3,2,5,3,5,4,3,6,7,6,6,5,6,7,7,7]

tets_list = full_snappy_triang(Dest)

#show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

print("Dest seq (index 8) is")
print(Dest)
print(' ')
print("Regular triangulation is")
show_triangulation(orb.Tetrahedra)
print(' ')
print("Isometry group of regular triangulation is")
print(orb.isometries())
print(' ')

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

#print((tet0.tilt(V1) + tet1.tilt(V2)).evaluate() > 0)

next_list = two_to_three(tets_list,tet0,F1)

#show_triangulation(next_list)

tet0 = CuspedOrbifold(next_list).Tetrahedra[0]

#print(tet0.tilt(V1).evaluate())
#print(tet0.tilt(V2).evaluate())

#print(check_2_to_3_possible([tet0],tet0,F1))

next_list = two_to_three([tet0],tet0,F1)

orb = CuspedOrbifold(next_list)

print("Canonical triangulation is")

show_triangulation(orb.Tetrahedra)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

#print(orb.is_canonical)
# It's canonical.

print(' ')
print("Isometry group of canonical triangulation is")
print(orb.isometries())
print(' ')
"""







"""
Dest = [0,1,2,3,2,2,0,2,1,0,1,1,4,3,3,0,3,4,4,4]

tets_list = full_snappy_triang(Dest)

show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

print(tet0.tilt(V0))
print(tet0.tilt(V1))
print(tet0.tilt(V2))
print(tet0.tilt(V3))
print(tet1.tilt(V2))
# Tilts of the faces of tet0 which are glued to themselves are negative. But tet0.tilt(V3) + tet1.tilt(V2) = 0,
# similarly tet1.tilt(Vi) + tet0.tilt(V3) = 0 for the other i's. This means we should attach a copy of tet0 to
# each face of tet1 and remove the connecting faces, giving a cube. The canonical decomposition is then this cube
# with a bunch of symmetries (all the ones preserving the inner tetrahedron) and the faces of the cube glued to
# themselves, carried over from tet0.
"""





"""
# One of the degree 14 covers in Fact 5.5.
Dest = [0,1,2,3, 2,4,0,2, 1,0,5,1, 6,3,3,0, 5,7,1,8, 4,2,7,9, 3,6,6,6, 7,5,4,10, 
11,9,11,4, 12,13,8,5, 13,12,13,7, 8,8,12,12, 9,11,10,11, 10,10,9,13] 

tets_list = full_snappy_triang(Dest)

show_triangulation(tets_list)
# It's a pretty interesting triangulation, definitely the most complicated so far.

orb = CuspedOrbifold(tets_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]

print(' ')
print(tet0.tilt(V0))
print(tet0.tilt(V1))
print(tet0.tilt(V2))
print(tet0.tilt(V3) + tet1.tilt(V2))
# We see that this is positive. Let's check the other tilts before anything else though.
print(' ')
print(tet1.tilt(V0) + tet2.tilt(V2))
print(tet2.tilt(V1))
print(' ')

# Alright, we're going to do a 2-3 move through face 3 of tet0 and face 2 of tet1. First
# check that it's possible.

print(check_2_to_3_possible(orb.Tetrahedra,tet0,F3))
print(' ')

# Returns true, so let's do it.

next_list = two_to_three(orb.Tetrahedra,tet0,F3)

show_triangulation(next_list)

orb = CuspedOrbifold(next_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

# From the gluing data, we know which tilts to check.

print(' ')
print((tet0.tilt(V0) + tet0.tilt(V3)).evaluate() < 0)
print(tet0.tilt(V1).evaluate() < 0)
print((tet0.tilt(V2) + tet1.tilt(V2)).evaluate() < 0)
print(tet0.tilt(V2) + tet1.tilt(V2))
print(tet1.tilt(V1))
print(tet1.tilt(V3))
print(tet1.tilt(V0))
print(' ')

print(check_2_to_3_possible(orb.Tetrahedra,tet0,F1))
print(check_2_to_3_possible(orb.Tetrahedra,tet0,F2))

# We're stuck. It's not canonical, but no 2-3 moves are allowed. It's possible the tilts are wrong,
# will do testing.

# IMPORTANT. As pointed out by Jason (see his email), there's a different kind of move you can do here.
# Take two copies of tet0, attach to the two faces of tet1 it should be attached to. That collection of
# three tetrahedra is mapped to itself by the symmetry of tet1. Can divide that up into some tet pieces.
# I'm not going to deal with this right now. Will work on other examples, will be interesting to see if
# this situation occurs again.
"""

"""
Following is the triangulation you get after doing a new move on the immediately above triangulation,
where we got stuck. See write-ups for explanantion of this move. It seems to work.
"""
"""
Dest = [0,1,2,3, 2,4,0,2, 1,0,5,1, 6,3,3,0, 5,7,1,8, 4,2,7,9, 3,6,6,6, 7,5,4,10, 
11,9,11,4, 12,13,8,5, 13,12,13,7, 8,8,12,12, 9,11,10,11, 10,10,9,13]

orb = dest_to_orb(Dest)
print("Dest seq (index 14) is")
print(Dest)
print(' ')
print("Regular triangulation is")
show_triangulation(orb.Tetrahedra)
print(' ')
print("Isometry group of regular triangulation is")
print(orb.isometries())
print(' ')


tet0 = Tetrahedron()
tet1 = Tetrahedron()
tet2 = Tetrahedron()

tet0.Index = 0
tet1.Index = 1
tet2.Index = 2

a = SquareRootCombination([(1,2)])
b = SquareRootCombination([(3,2)])
tet0.fill_edge_params(ComplexSquareRootCombination(a,b))

a = SquareRootCombination.Zero()
b = SquareRootCombination([(3,Fraction(1,6))])
tet1.fill_edge_params(ComplexSquareRootCombination(a,b))

a = SquareRootCombination([(1,Fraction(21,18))])
b = SquareRootCombination([(3,Fraction(1,6))])
tet2.fill_edge_params(ComplexSquareRootCombination(a,b))

tet0.attach(F0,tet1,(1,0,2,3))
tet0.attach(F2,tet1,(1,0,2,3))
tet1.attach(F0,tet2,(3,1,2,0))
tet1.attach(F3,tet2,(3,1,2,0))

tet0.Symmetries = [Perm4((0,1,2,3)),Perm4((1,0,3,2))]
tet2.Symmetries = [Perm4((0,1,2,3)),Perm4((2,3,0,1))]

orb = CuspedOrbifold([tet0,tet1,tet2])

print("Canonical triangultion is")
show_triangulation(orb.Tetrahedra)
print(' ')
print("Isometry group of canonical triangulation is")
print(orb.isometries())



print(orb.is_canonical)

print(tet0.tilt(V0) + tet1.tilt(V1))
print(tet0.tilt(V2) + tet1.tilt(V2))
print(tet1.tilt(V0) + tet2.tilt(V3))
print(tet1.tilt(V3) + tet2.tilt(V0))

# Yay! It seems like this is canonical. So we had to do a crazy new move, but it worked. 
"""




"""
# Another degree 14 cover.
Dest = [1,1,1,2, 0,0,0,3, 4,4,5,0, 6,7,6,1, 2,8,2,6, 8,2,9,7, 3,3,10,4, 10,11,3,5, 5,12,4,11, 
12,5,12,12, 7,6,13,13, 13,13,7,8, 9,9,8,9, 11,10,11,10]

tets_list = full_snappy_triang(Dest)

show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]
tet3 = orb.Tetrahedra[3]

print(tet0.tilt(V2) + tet1.tilt(V2))
print(tet0.tilt(V0) + tet2.tilt(V2))
print(tet2.tilt(V0) + tet3.tilt(V2))

# All negative, so it's canonical. That's good, because I don't think any 2-3 moves were possible here.
"""



"""
# The one triangulation with 9 tetrahedra.

Dest = [0,0,0,1, 2,3,4,0, 1,5,6,2, 6,6,1,4, 5,1,5,3, 4,4,2,7, 3,2,3,8, 8,8,8,5, 7,7,7,6]

tets_list = full_snappy_triang(Dest)

show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

print(tet0.tilt(V0) + tet1.tilt(V2))
print(tet0.tilt(V2))

# Both are negative, so this is canonical.
"""


"""
# A triangulation with 12 tetrahedra.
Dest = [0,1,1,2, 1,0,0,3, 4,5,5,0, 6,7,7,1, 2,8,8,4, 8,2,2,7,
 3,9,9,6, 9,3,3,5, 5,4,4,10, 7,6,6,11, 11,10,10,8, 10,11,11,9]

orb = dest_to_orb(Dest)

show_triangulation(orb.Tetrahedra)
print(orb.is_canonical)
print(orb.DestSeq)

tet0 = orb.Tetrahedra[0]

print(tet0.tilt(V0) + tet0.tilt(V1))
print((tet0.tilt(V0) + tet0.tilt(V1)).evaluate())
print(tet0.tilt(V2))
print(tet0.tilt(V3))

print(orb.is_canonical)

# It's canonical.
"""


"""
# Next with 12 tets.
Dest = [0,1,2,0, 2,3,0,4, 1,0,5,6, 5,7,1,5, 8,6,6,1, 3,2,7,3,
 9,4,4,2, 7,5,3,10, 4,9,9,9, 6,8,8,8, 11,10,10,7, 10,11,11,11]

tets_list = full_snappy_triang(Dest)

show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]

print(tet0.tilt(V2))
# this is positive

print(tet0.tilt(V0) + tet1.tilt(V2))

print(tet1.tilt(V3) + tet2.tilt(V2))
# this is also positive

print(tet2.tilt(V0))

# We should be able to do a 2-3 move through either of the two bad faces.

print(check_2_to_3_possible(orb.Tetrahedra,tet0,F2))
print(check_2_to_3_possible(orb.Tetrahedra,tet1,F3))

next_list = two_to_three(orb.Tetrahedra,tet0,F2)

show_triangulation(next_list)

orb = CuspedOrbifold(next_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]

print(tet0.tilt(V0) + tet0.tilt(V3))
print(tet0.tilt(V1) + tet1.tilt(V2))
print(tet1.tilt(V3) + tet2.tilt(V2))
# This is still a problem.

print(check_2_to_3_possible(orb.Tetrahedra,tet1,F3))

next_list = two_to_three(orb.Tetrahedra,tet1,F3)

show_triangulation(next_list)

orb = CuspedOrbifold(next_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

print(tet0.tilt(V2))
print(tet0.tilt(V1) + tet1.tilt(V1))
# These are both positive, but we can't do a 2-3 to F1 of tet0, so let's do it to F2 of tet0.

next_list = two_to_three(orb.Tetrahedra,tet0,F2)

show_triangulation(next_list)

orb = CuspedOrbifold(next_list)
tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]

print(tet0.tilt(V0) + tet1.tilt(V1))
# This is still positive, but we're no longer able to do any 2-3 moves because of the symmetry
# of tet0. That means we have to do a different kind of move, as in a previous example. For now we leave
# this example alone. 
"""



"""
# WARNING. There's a problem with the following computations. Did 2-3 moves which were illegal
# because of large dihedral angles (I think). Forgot to check beforehand. So we didn't succeed in getting to the canonical decomp,
# in fact I don't think we can with only 2-3 moves, in this situation.
# Let's do the same one as above, except do the other 2-3 move option at the first step
Dest = [0,1,2,0, 2,3,0,4, 1,0,5,6, 5,7,1,5, 8,6,6,1, 3,2,7,3,
 9,4,4,2, 7,5,3,10, 4,9,9,9, 6,8,8,8, 11,10,10,7, 10,11,11,11]

orb = dest_to_orb(Dest)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]

show_triangulation(orb.Tetrahedra)

next_orb = CuspedOrbifold(two_to_three(orb.Tetrahedra,tet1,F3))

tet0 = next_orb.Tetrahedra[0]
tet1 = next_orb.Tetrahedra[1]

show_triangulation(next_orb.Tetrahedra)


print(next_orb.is_canonical)
# returns false. let's check where it's bad.

print(tet0.tilt(V2))
# this is positive

print(tet0.tilt(V1) + tet1.tilt(V0))
# this is zero

print(tet1.tilt(V2))
# this is positive

print(tet0.tilt(V0) + tet0.tilt(V3))

next_orb = CuspedOrbifold(two_to_three(next_orb.Tetrahedra,tet0,F2))

show_triangulation(next_orb.Tetrahedra)

tet0 = next_orb.Tetrahedra[0]
tet1 = next_orb.Tetrahedra[1]
tet2 = next_orb.Tetrahedra[2]

print(next_orb.is_canonical)
# returns false still

print(tet0.tilt(V0) + tet1.tilt(V1))
# positive

print(tet0.tilt(V1) + tet2.tilt(V0))

print(tet1.tilt(V3))

print(tet2.tilt(V2))
# positive

# Can't do a 2-3 move through F0 of tet0 like in the previous examples, could do the "new" move though.
# But we can do a 2-3 move though F2 of tet2. Let's do that.

next_orb = CuspedOrbifold(two_to_three(next_orb.Tetrahedra,tet2,F2))

show_triangulation(next_orb.Tetrahedra)

print(next_orb.is_canonical)
# apparently this is canonical.
"""

"""
# the 8-tet dest seq whose isometry group jason wants to know.
Dest = [1,2,2,1,0,3,3,0,3,0,0,4,2,1,1,5,6,5,5,2,7,4,4,3,4,7,7,7,5,6,6,6]

orb = dest_to_orb(Dest)

print("The data of the triangulation is")

print(Dest)

show_triangulation(orb.Tetrahedra)

print("The isometries are")

print(orb.isometries())
"""


"""
Now let's do some testing for my canonize function, found in The_Algorithm.py.
"""
"""
Dest = [0,1,1,0,1,0,0,2,3,2,2,1,2,3,3,3]
orb = dest_to_orb(Dest)
canonical_orb = canonize(orb)
print(canonical_orb.DestSeq)
print(canonical_orb.PachnerPath)
show_triangulation(canonical_orb.Tetrahedra)
"""
"""
Dest = [0,0,0,0]
orb = dest_to_orb(Dest)
canonical_orb = canonize(orb)
"""
"""
Dest = [0,1,2,3,2,2,0,2,1,0,1,1,4,3,3,0,3,4,4,4]
orb = dest_to_orb(Dest)
canonical_orb = canonize(orb)
"""
"""
Dest = [0,1,2,1,2,3,0,0,1,0,4,2,4,5,1,4,3,2,5,3,5,4,3,6,7,6,6,5,6,7,7,7]
orb = dest_to_orb(Dest)
canonical_orb = canonize(orb)
print(canonical_orb.DestSeq)
print(canonical_orb.PachnerPath)
show_triangulation(canonical_orb.Tetrahedra)
"""
"""
Dest = [0,1,2,3, 2,4,0,2, 1,0,5,1, 6,3,3,0, 5,7,1,8, 4,2,7,9, 3,6,6,6, 7,5,4,10, 
11,9,11,4, 12,13,8,5, 13,12,13,7, 8,8,12,12, 9,11,10,11, 10,10,9,13]
orb = dest_to_orb(Dest)
canonical_orb = canonize(orb)
# this is the one where you get stuck if you only do 2-3 moves. And canonize gives the
# correct error message.
"""
"""
Dest = [0,1,2,0, 2,3,0,4, 1,0,5,6, 5,7,1,5, 8,6,6,1, 3,2,7,3,
 9,4,4,2, 7,5,3,10, 4,9,9,9, 6,8,8,8, 11,10,10,7, 10,11,11,11]
orb = dest_to_orb(Dest)
canonical_orb = canonize(orb)
print(canonical_orb.DestSeq)
print(canonical_orb.PachnerPath)
show_triangulation(canonical_orb.Tetrahedra)
"""


"""
Dest = [0,1,2,0, 2,3,0,4, 1,0,5,6, 5,7,1,5, 8,6,6,1, 3,2,7,3,
 9,4,4,2, 7,5,3,10, 4,9,9,9, 6,8,8,8, 11,10,10,7, 10,11,11,11]
orb = dest_to_orb(Dest)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]

print(tet0.horotriangles)
print(tet1.horotriangles)
print(tet2.horotriangles)

next_orb = CuspedOrbifold(two_to_three(orb.Tetrahedra,tet0,F2))

tet0 = next_orb.Tetrahedra[0]
tet1 = next_orb.Tetrahedra[1]
tet2 = next_orb.Tetrahedra[2]

print(tet0.horotriangles)
print(tet1.horotriangles)
print(tet2.horotriangles)
"""

"""
with open("newly_canonical_dest_seqs.json", "r") as read_file:
    newly_canonical_dest_seqs = json.load(read_file)

for dest in newly_canonical_dest_seqs:
    if len(dest) == 32:
        print(dest)
"""

"""
with open("stuck_dest_seqs.json", "r") as read_file:
    stuck_dest_seqs = json.load(read_file)

Dest = stuck_dest_seqs[5]
orb = dest_to_orb(Dest)
show_triangulation(orb.Tetrahedra)
tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]

orb.arrow_two_to_three(F3,tet0)
show_triangulation(orb.Tetrahedra)
tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
print((tet0.tilt(V2)+tet1.tilt(V1)).evaluate() > 0)
print(check_2_to_3_possible(orb.Tetrahedra,tet0,F2))
#this is where we get stuck, want to do a 2-3 move but there's a flat quad.
"""

a = 1
b = 2
if (a != 1
    or b == 2):
    print('ye')

"""
For dest = stuck_dest_seqs[5], we can do one 2-3 move then we get stuck, can't do more 2-3 moves
because of large dihedral angles. But one angle sum is actually pi, resulting in a flat quad face,
and it turns out we can do a modified 2-3 move, modified in the sense that we can just remove the flat
tetrahedron resulting from the 2-3 move. See my notebook for more info. I constructed by hand the resulting 
triangulation, here it is.
"""
"""
tet0 = Tetrahedron()
tet0.Index = 0
tet1 = Tetrahedron()
tet1.Index = 1

tet0.attach(F0,tet1,(1,0,2,3))
tet0.attach(F1,tet1,(1,3,0,2))
tet0.attach(F2,tet1,(1,3,0,2))
tet0.attach(F3,tet1,(1,3,0,2))

a = SquareRootCombination.One()
b = SquareRootCombination([(3,Fraction(1,3))])
tet0.fill_edge_params(ComplexSquareRootCombination(a,b))

a = SquareRootCombination.Zero()
b = SquareRootCombination([(3,Fraction(1,3))])
tet1.fill_edge_params(ComplexSquareRootCombination(a,b))

orb = CuspedOrbifold([tet0,tet1])

print(orb.is_canonical)
show_triangulation(orb.Tetrahedra)

print(orb.Vertices)
print(orb.Vertices[0].Corners)
print(orb.Vertices[1].Corners)
"""

"""
tet0 = Tetrahedron()
tet0.Index = 0
tet1 = Tetrahedron()
tet1.Index = 1

tet0.attach(F0,tet1,(1,3,0,2))
tet0.attach(F1,tet1,(1,3,0,2))
tet0.attach(F2,tet1,(1,3,0,2))
tet0.attach(F3,tet1,(1,3,0,2))

a = SquareRootCombination.One()
b = SquareRootCombination([(3,Fraction(1,3))])
tet0.fill_edge_params(ComplexSquareRootCombination(a,b))

a = SquareRootCombination.Zero()
b = SquareRootCombination([(3,Fraction(1,3))])
tet1.fill_edge_params(ComplexSquareRootCombination(a,b))

orb = CuspedOrbifold([tet0,tet1])

print(orb.is_canonical)
show_triangulation(orb.Tetrahedra)

print(orb.Vertices)
print(orb.Vertices[0].Corners)
print(orb.Vertices[1].Corners)
print(orb.Vertices[2].Corners)
print(orb.Vertices[3].Corners)
"""






"""
For each dest seq in Dests, tries to get the isometry group of the canonical decomp and prints.
At this point, can only get the canonical decomp using 2-3 moves. These are dest seqs which we know
cover the 4-tet dest seq, which are important because they might cover stuff which is not in C_main.
"""
"""
Dests = [[1, 2, 2, 1, 0, 3, 3, 0, 3, 0, 0, 4, 2, 1, 1, 5, 6, 5, 5, 2, 7, 4, 4, 3, 4, 7, 7, 7, 5, 6, 6, 6],[0, 1, 1, 2, 1, 0, 0, 3, 4, 5, 5, 0, 6, 7, 7, 1, 2, 8, 8, 4, 8, 2, 2, 7, 3, 9, 9, 6, 9, 3, 3, 5, 5, 4, 4, 10, 7, 6, 6, 11, 11, 10, 10, 8, 10, 11, 11, 9]
,[0, 1, 2, 0, 2, 3, 0, 4, 1, 0, 5, 6, 5, 7, 1, 5, 8, 6, 6, 1, 3, 2, 7, 3, 9, 4, 4, 2, 7, 5, 3, 10, 4, 9, 9, 9, 6, 8, 8, 8, 11, 10, 10, 7, 10, 11, 11, 11]
,[0, 1, 2, 3, 2, 3, 0, 4, 1, 0, 5, 6, 5, 7, 1, 0, 8, 4, 4, 1, 3, 2, 7, 5, 9, 10, 10, 2, 7, 5, 3, 10, 4, 8, 8, 9, 6, 11, 11, 8, 11, 6, 6, 7, 10, 9, 9, 11]
,[0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 3, 6, 4, 2, 7, 0, 3, 7, 1, 4, 8, 9, 9, 1, 10, 6, 6, 2, 7, 3, 4, 9, 5, 11, 11, 10, 11, 5, 5, 7, 6, 10, 10, 8, 9, 8, 8, 11]
,[0, 1, 2, 3, 2, 4, 0, 2, 1, 0, 5, 1, 6, 3, 3, 0, 5, 7, 1, 8, 4, 2, 7, 9, 3, 6, 6, 6, 7, 5, 4, 10, 11, 9, 12, 4, 13, 12, 8, 5, 14, 15, 15, 7, 8, 16, 13, 13, 16, 8, 9, 15, 9, 11, 16, 11, 10, 17, 17, 14, 17, 10, 10, 12, 12, 13, 11, 18, 15, 14, 14, 19, 19, 18, 18, 16, 18, 19, 19, 17]
,[0, 1, 2, 3, 2, 4, 0, 1, 1, 0, 5, 6, 7, 8, 9, 0, 5, 10, 1, 9, 4, 2, 10, 11, 6, 12, 13, 2, 3, 14, 15, 7, 15, 9, 3, 13, 14, 3, 8, 4, 10, 5, 4, 16, 17, 18, 19, 5, 13, 20, 6, 19, 12, 6, 16, 8, 9, 15, 7, 17, 8, 7, 14, 21, 20, 13, 22, 10, 11, 21, 23, 14, 23, 19, 11, 22, 21, 11, 18, 12, 16, 22, 12, 20, 19, 23, 17, 15, 22, 16, 20, 18, 18, 17, 21, 23]
,[1, 2, 2, 2, 0, 3, 3, 4, 3, 0, 0, 0, 2, 1, 1, 5, 5, 6, 7, 1, 4, 8, 9, 3, 9, 10, 4, 9, 8, 4, 11, 8, 7, 12, 5, 7, 6, 5, 13, 6, 13, 14, 6, 15, 12, 7, 14, 16, 11, 17, 8, 18, 10, 9, 17, 19, 17, 11, 10, 17, 20, 19, 19, 10, 21, 18, 18, 11, 14, 13, 12, 14, 22, 16, 16, 12, 23, 15, 15, 13, 15, 23, 23, 23, 16, 22, 22, 22, 18, 21, 21, 21, 19, 20, 20, 20]
,[1, 2, 2, 3, 0, 4, 4, 5, 4, 0, 0, 6, 7, 8, 8, 0, 2, 1, 1, 9, 10, 11, 11, 1, 12, 13, 13, 2, 3, 14, 14, 10, 14, 3, 3, 13, 15, 16, 16, 4, 5, 17, 17, 7, 17, 5, 5, 16, 6, 18, 18, 15, 18, 6, 6, 8, 8, 7, 7, 19, 9, 20, 20, 12, 20, 9, 9, 11, 11, 10, 10, 21, 13, 12, 12, 22, 22, 21, 21, 14, 16, 15, 15, 23, 23, 19, 19, 17, 19, 23, 23, 18, 21, 22, 22, 20]
,[1, 2, 2, 3, 0, 4, 4, 5, 4, 0, 0, 6, 7, 8, 8, 0, 2, 1, 1, 9, 10, 11, 12, 1, 13, 14, 14, 2, 3, 15, 15, 10, 15, 3, 3, 14, 16, 17, 18, 4, 5, 18, 19, 7, 19, 20, 5, 18, 18, 5, 16, 17, 6, 21, 21, 16, 21, 6, 6, 8, 8, 7, 7, 20, 9, 12, 22, 13, 22, 23, 9, 12, 12, 9, 10, 11, 11, 10, 23, 22, 23, 22, 11, 15, 14, 13, 13, 23, 17, 16, 20, 19, 20, 19, 17, 21]
,[1, 2, 2, 3, 0, 4, 4, 5, 4, 0, 0, 6, 7, 8, 8, 0, 2, 1, 1, 9, 10, 11, 12, 1, 13, 14, 14, 2, 3, 15, 15, 10, 15, 3, 3, 14, 16, 17, 18, 4, 5, 19, 17, 7, 17, 16, 5, 18, 19, 5, 20, 17, 6, 21, 21, 16, 21, 6, 6, 8, 8, 7, 7, 20, 9, 22, 11, 13, 11, 10, 9, 12, 22, 9, 23, 11, 12, 23, 10, 22, 23, 12, 22, 15, 14, 13, 13, 23, 18, 20, 16, 19, 20, 18, 19, 21]
,[0, 1, 2, 3, 2, 4, 0, 1, 1, 0, 5, 6, 7, 8, 9, 0, 5, 10, 1, 9, 4, 2, 10, 11, 6, 12, 13, 2, 3, 14, 15, 7, 15, 16, 3, 13, 14, 3, 17, 4, 10, 5, 4, 18, 19, 17, 16, 5, 13, 20, 6, 16, 12, 6, 21, 8, 9, 22, 7, 19, 8, 7, 23, 23, 23, 11, 8, 12, 22, 9, 11, 24, 25, 24, 24, 10, 11, 23, 22, 14, 21, 26, 12, 21, 20, 13, 26, 20, 17, 19, 14, 27, 16, 15, 19, 15, 28, 18, 18, 17, 18, 28, 28, 25, 26, 21, 20, 29, 30, 27, 27, 22, 24, 25, 25, 30, 31, 29, 29, 26, 27, 30, 30, 28, 29, 31, 31, 31]
,[0, 1, 2, 3, 2, 4, 0, 2, 1, 0, 5, 1, 6, 3, 3, 0, 5, 7, 1, 8, 4, 2, 7, 9, 3, 6, 6, 6, 7, 5, 4, 10, 11, 9, 12, 4, 13, 14, 8, 5, 15, 16, 17, 7, 8, 18, 13, 13, 18, 8, 19, 16, 9, 11, 20, 11, 20, 21, 9, 17, 10, 17, 22, 15, 22, 23, 10, 12, 17, 10, 15, 14, 12, 24, 11, 25, 24, 12, 21, 26, 14, 13, 25, 20, 25, 19, 14, 22, 16, 15, 23, 21, 23, 22, 16, 27, 19, 25, 18, 28, 21, 20, 24, 18, 29, 27, 27, 19, 30, 26, 26, 23, 31, 28, 28, 24, 26, 30, 30, 31, 27, 29, 29, 30, 28, 31, 31, 29]
,[0, 1, 1, 2, 1, 0, 0, 3, 4, 5, 6, 0, 7, 8, 9, 1, 2, 10, 11, 4, 11, 12, 2, 9, 10, 2, 13, 8, 3, 14, 15, 7, 15, 16, 3, 6, 14, 3, 17, 5, 6, 18, 4, 19, 5, 4, 20, 21, 20, 22, 5, 23, 18, 6, 22, 24, 9, 25, 7, 26, 8, 7, 27, 28, 27, 29, 8, 30, 25, 9, 29, 31, 13, 32, 10, 18, 28, 21, 33, 10, 12, 11, 32, 20, 26, 33, 19, 11, 32, 13, 12, 29, 23, 31, 31, 12, 24, 30, 30, 13, 17, 34, 14, 25, 21, 28, 35, 14, 16, 15, 34, 27, 19, 35, 26, 15, 34, 17, 16, 22, 30, 24, 24, 16, 31, 23, 23, 17, 22, 20, 18, 33, 35, 19, 21, 32, 29, 27, 25, 35, 33, 26, 28, 34]
,[0, 1, 1, 2, 1, 0, 0, 3, 4, 5, 6, 0, 7, 8, 9, 1, 2, 10, 11, 4, 11, 12, 2, 9, 10, 2, 13, 8, 3, 14, 15, 7, 15, 16, 3, 6, 14, 3, 17, 5, 6, 18, 4, 19, 5, 4, 20, 21, 20, 22, 5, 23, 18, 6, 22, 24, 9, 25, 7, 26, 8, 7, 27, 28, 27, 29, 8, 30, 25, 9, 29, 31, 13, 32, 10, 20, 28, 21, 21, 10, 12, 11, 32, 18, 26, 19, 19, 11, 32, 13, 12, 29, 24, 31, 31, 12, 23, 30, 30, 13, 17, 33, 14, 27, 21, 28, 28, 14, 16, 15, 33, 25, 19, 26, 26, 15, 33, 17, 16, 22, 31, 24, 24, 16, 30, 23, 23, 17, 22, 20, 18, 34, 29, 27, 25, 35, 35, 34, 34, 32, 34, 35, 35, 33]
,[0, 1, 2, 0, 2, 3, 0, 4, 1, 0, 5, 6, 5, 7, 1, 8, 9, 6, 10, 1, 3, 2, 7, 11, 12, 13, 4, 2, 7, 5, 3, 14, 15, 16, 17, 3, 4, 18, 12, 12, 18, 4, 19, 16, 20, 17, 21, 5, 6, 9, 22, 9, 22, 23, 6, 21, 24, 25, 25, 7, 8, 26, 27, 20, 27, 28, 8, 10, 26, 8, 11, 25, 10, 29, 9, 29, 29, 10, 23, 27, 11, 30, 26, 15, 30, 11, 28, 13, 13, 12, 31, 31, 31, 19, 13, 30, 14, 32, 32, 24, 32, 14, 14, 17, 17, 20, 15, 33, 16, 15, 34, 19, 34, 21, 16, 34, 19, 31, 18, 18, 21, 34, 20, 23, 23, 22, 29, 22, 25, 24, 24, 35, 35, 33, 33, 26, 28, 27, 30, 28, 33, 35, 35, 32]
,[0, 1, 2, 0, 2, 3, 0, 4, 1, 0, 5, 6, 5, 7, 1, 8, 9, 6, 10, 1, 3, 2, 7, 11, 12, 13, 4, 2, 7, 5, 3, 14, 15, 16, 17, 3, 4, 18, 12, 12, 18, 4, 19, 16, 20, 21, 22, 5, 6, 9, 23, 9, 23, 24, 6, 22, 25, 26, 27, 7, 8, 17, 28, 20, 28, 29, 8, 10, 17, 8, 15, 26, 10, 30, 9, 30, 30, 10, 24, 28, 11, 31, 21, 15, 21, 20, 11, 27, 31, 11, 32, 13, 13, 12, 33, 33, 33, 19, 13, 31, 14, 34, 35, 25, 35, 27, 14, 17, 34, 14, 26, 21, 16, 15, 29, 19, 29, 28, 16, 29, 19, 33, 18, 18, 22, 32, 20, 24, 32, 22, 31, 32, 24, 23, 30, 23, 27, 35, 25, 34, 26, 25, 34, 35]
,[0, 1, 2, 3, 2, 4, 0, 1, 1, 0, 5, 6, 7, 8, 9, 0, 5, 10, 1, 9, 4, 2, 10, 11, 6, 12, 13, 2, 3, 14, 15, 7, 15, 16, 3, 13, 14, 3, 17, 4, 10, 5, 4, 18, 19, 20, 21, 5, 13, 22, 6, 21, 12, 6, 23, 8, 9, 24, 7, 19, 8, 7, 25, 26, 25, 27, 8, 28, 24, 9, 27, 29, 30, 29, 31, 10, 11, 26, 32, 14, 32, 21, 11, 31, 26, 11, 20, 12, 23, 28, 12, 33, 22, 13, 28, 23, 17, 34, 14, 24, 16, 15, 34, 25, 21, 32, 19, 15, 34, 17, 16, 35, 28, 23, 22, 16, 29, 30, 18, 17, 18, 35, 29, 30, 35, 18, 33, 20, 20, 19, 26, 34, 33, 31, 35, 22, 27, 25, 24, 32, 31, 33, 30, 27]
,[0, 1, 2, 3, 2, 4, 0, 2, 1, 0, 5, 1, 6, 3, 3, 0, 5, 7, 1, 8, 4, 2, 7, 9, 3, 6, 6, 6, 7, 5, 4, 10, 11, 9, 12, 4, 13, 14, 8, 5, 15, 16, 17, 7, 8, 18, 13, 13, 18, 8, 19, 16, 9, 11, 20, 11, 20, 21, 9, 17, 10, 22, 23, 15, 23, 24, 10, 12, 22, 10, 25, 14, 12, 26, 11, 26, 26, 12, 21, 23, 14, 13, 27, 27, 27, 19, 14, 22, 17, 28, 15, 21, 16, 15, 29, 19, 29, 30, 16, 29, 28, 17, 30, 28, 19, 27, 18, 18, 21, 20, 26, 20, 25, 31, 22, 25, 24, 23, 31, 24, 31, 25, 24, 32, 30, 29, 28, 33, 34, 33, 33, 30, 35, 32, 32, 31, 32, 35, 35, 35, 33, 34, 34, 34]
,[0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 6, 2, 7, 8, 9, 0, 6, 10, 1, 11, 5, 12, 13, 1, 4, 2, 10, 8, 3, 14, 15, 7, 15, 16, 3, 6, 14, 3, 17, 12, 10, 6, 4, 18, 19, 20, 21, 4, 13, 22, 5, 9, 12, 5, 23, 20, 9, 24, 7, 25, 8, 7, 26, 19, 26, 27, 8, 28, 24, 9, 27, 29, 30, 31, 28, 10, 11, 32, 25, 15, 25, 21, 11, 13, 32, 11, 20, 31, 23, 29, 12, 22, 22, 13, 29, 33, 17, 34, 14, 24, 20, 19, 32, 14, 16, 15, 34, 26, 34, 17, 16, 35, 28, 18, 30, 16, 29, 23, 22, 17, 18, 28, 35, 30, 35, 33, 18, 21, 21, 25, 19, 34, 33, 35, 31, 23, 27, 26, 24, 32, 31, 30, 33, 27]
,[0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 6, 7, 8, 9, 9, 0, 6, 10, 1, 11, 12, 13, 14, 1, 4, 2, 10, 15, 16, 17, 13, 2, 3, 18, 18, 8, 18, 3, 3, 13, 10, 6, 4, 19, 20, 21, 21, 4, 5, 22, 23, 16, 23, 7, 5, 9, 22, 5, 19, 21, 24, 25, 25, 6, 7, 23, 26, 12, 26, 19, 7, 25, 9, 8, 8, 27, 28, 14, 17, 10, 11, 29, 29, 24, 29, 11, 11, 14, 14, 28, 12, 30, 13, 12, 16, 31, 15, 32, 32, 20, 32, 15, 15, 17, 17, 16, 28, 33, 31, 27, 27, 18, 19, 26, 22, 28, 21, 20, 20, 34, 34, 33, 33, 22, 27, 31, 31, 23, 25, 24, 24, 35, 35, 30, 30, 26, 30, 35, 35, 29, 33, 34, 34, 32]
,[0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 6, 7, 8, 9, 9, 0, 6, 10, 1, 11, 12, 13, 14, 1, 4, 2, 10, 15, 16, 17, 13, 2, 3, 18, 18, 8, 18, 3, 3, 13, 10, 6, 4, 19, 20, 21, 22, 4, 5, 23, 24, 16, 24, 7, 5, 9, 23, 5, 25, 21, 26, 27, 28, 6, 7, 24, 29, 12, 29, 25, 7, 28, 9, 8, 8, 30, 25, 29, 23, 10, 11, 31, 32, 26, 32, 28, 11, 14, 31, 11, 27, 29, 14, 19, 12, 27, 13, 12, 16, 33, 19, 14, 17, 25, 15, 34, 35, 20, 35, 22, 15, 23, 34, 15, 21, 17, 17, 16, 19, 22, 33, 30, 30, 18, 22, 35, 20, 34, 21, 20, 34, 35, 30, 33, 33, 24, 28, 32, 26, 31, 27, 26, 31, 32]
,[0, 1, 2, 3, 2, 4, 0, 5, 1, 0, 6, 7, 8, 9, 10, 0, 6, 11, 1, 12, 13, 14, 15, 1, 4, 2, 11, 16, 17, 15, 18, 2, 3, 19, 20, 8, 20, 12, 3, 18, 19, 3, 16, 14, 11, 6, 4, 21, 22, 23, 9, 4, 5, 24, 25, 17, 25, 21, 5, 10, 24, 5, 7, 23, 26, 10, 23, 6, 7, 27, 24, 13, 27, 7, 21, 9, 10, 26, 8, 28, 9, 8, 22, 29, 30, 18, 14, 11, 12, 20, 31, 26, 31, 16, 12, 15, 15, 17, 13, 32, 14, 13, 30, 33, 16, 31, 19, 22, 18, 30, 17, 34, 33, 29, 29, 19, 34, 28, 28, 20, 21, 25, 27, 30, 23, 22, 26, 35, 35, 32, 32, 24, 28, 34, 34, 25, 29, 33, 33, 27, 32, 35, 35, 31]]


for i in range(len(Dests)):
    print(Dests[i])
    orb = dest_to_orb(Dests[i])
    if orb.is_canonical:
        print(' ')
        print("Orb with regular tet triangulation is already canonical. Its triangulation is")
        print(' ')
        show_triangulation(orb.Tetrahedra)
        print(' ')
        print("Its isometry group is")
        print(" ")
        print(orb.isometries())
    else:
        print(' ')
        print("Orb with regular tet triangulation is not canonical. The number of symmetries of the regular triangulation is")
        print(' ')
        print(len(orb.isometries()))
        print(' ')
        orb = canonize(orb)
        if orb != None:
            print("Orb can be canonized with 2-3 moves.")
            TiltSumZero = False
            for tet1 in orb.Tetrahedra:
                for face1 in TwoSubsimplices:
                    if tet1.Neighbor[face1] != None:
                        tet2 = tet1.Neighbor[face1]
                        face2 = tet1.Gluing[face1].image(face1)
                        if (tet1.tilt(comp(face2)) + tet2.tilt(comp(face2))).evaluate() == 0:
                            TiltSumZero = True
            if TiltSumZero:
                print("Canonical decomp has non-tetrahedral cells. Will not compute isometry group.")
            else:
                print("Canonical triangulation is")
                print(' ')
                show_triangulation(orb.Tetrahedra)
                print(' ')
                print("Isometry group is")
                print(" ")
                print(orb.isometries())
    print(" ")
    print(" ")
    print("------------------------------------------------------")
    print(' ')
    print(' ')
"""
"""
Dest = [0,1,1,0,1,0,0,2,3,2,2,1,2,3,3,3]
orb = dest_to_orb(Dest)
show_triangulation(orb.Tetrahedra)
tet0 = orb.Tetrahedra[0]

orb.arrow_two_to_three(F2,tet0)
show_triangulation(orb.Tetrahedra)

tet0 = orb.Tetrahedra[0]
orb.three_to_two(tet0.Class[E01])
show_triangulation(orb.Tetrahedra)
#passes this test
"""
"""
Dest = [0,1,2,1,2,3,0,0,1,0,4,2,4,5,1,4,3,2,5,3,5,4,3,6,7,6,6,5,6,7,7,7]
orb = dest_to_orb(Dest)
show_triangulation(orb.Tetrahedra)
tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

orb.arrow_two_to_three(F1,tet0)
show_triangulation(orb.Tetrahedra)
tet0 = orb.Tetrahedra[0]

orb.three_to_two(tet0.Class[E01])
show_triangulation(orb.Tetrahedra)
#Also works for this
"""

