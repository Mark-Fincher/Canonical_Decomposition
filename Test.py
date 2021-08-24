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


"""
Dest = [0,0,0,0]

tets_list = full_snappy_triang(Dest)

show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

tet0 = orb.Tetrahedra[0]

print(tet0.horotriangles[V0].area)

print(tet0.tilt(V0))
print(tet0.tilt(V1))
print(tet0.tilt(V2))
print(tet0.tilt(V3))
"""



"""
Dest = [0,1,1,0,1,0,0,2,3,2,2,1,2,3,3,3]

tets_list = full_snappy_triang(Dest)

#show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

tet0 = orb.Tetrahedra[0]

print(tet0.tilt(V0))
print(tet0.tilt(V1))
print(tet0.tilt(V2))
print(tet0.tilt(V3))

# tet0.tilt(V2) is positive, and F2 is glued to itself, so that face is bad. Do 2-3 move through it.

#print(check_2_to_3_possible(orb.Tetrahedra,tet0,F2))

next_list = two_to_three(orb.Tetrahedra,orb.Tetrahedra[0],F2)

show_triangulation(next_list)

orb = CuspedOrbifold(next_list)

tet0 = orb.Tetrahedra[0]

print(tet0.horotriangles[V0].lengths)
print(tet0.horotriangles[V0].circumradius)
print(tet0.horotriangles[V1].lengths)
print(tet0.horotriangles[V1].circumradius)
print(tet0.horotriangles[V2].lengths)
print(tet0.horotriangles[V2].circumradius)
print(tet0.horotriangles[V3].lengths)
print(tet0.horotriangles[V3].circumradius)

print(tet0.tilt(V0))
print(tet0.tilt(V1))
print(tet0.tilt(V2))
print(tet0.tilt(V3))

# from checking tilt sums, this is canonical.
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
Dest = [0,1,2,1,2,3,0,0,1,0,4,2,4,5,1,4,3,2,5,3,5,4,3,6,7,6,6,5,6,7,7,7]

tets_list = full_snappy_triang(Dest)

show_triangulation(tets_list)

orb = CuspedOrbifold(tets_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

print((tet0.tilt(V1) + tet1.tilt(V2)).evaluate() > 0)

next_list = two_to_three(tets_list,tet0,F1)

show_triangulation(next_list)

tet0 = CuspedOrbifold(next_list).Tetrahedra[0]

print(tet0.tilt(V1).evaluate())
print(tet0.tilt(V2).evaluate())

print(check_2_to_3_possible([tet0],tet0,F1))

next_list = two_to_three([tet0],tet0,F1)

show_triangulation(next_list)

orb = CuspedOrbifold(next_list)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]

print(tet0.tilt(V0) + tet1.tilt(V1))
print((tet0.tilt(V0) + tet1.tilt(V1)).evaluate() > 0)

# returns false, so that face is fine. Now check the others.

print(tet0.tilt(V1) + tet1.tilt(V2))
print((tet0.tilt(V1) + tet1.tilt(V2)).evaluate() > 0)

print(tet1.tilt(V0))
print(tet1.tilt(V0).evaluate() > 0)

print(tet1.tilt(V3))
print(tet1.tilt(V3).evaluate() > 0)

# all return false, so this is canonical!
"""




"""
# One of the degree 14 covers.
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

show_triangulation(orb.Tetrahedra)
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
# Let's do the same one as above, except do the other 2-3 move option at the first step
Dest = [0,1,2,0, 2,3,0,4, 1,0,5,6, 5,7,1,5, 8,6,6,1, 3,2,7,3,
 9,4,4,2, 7,5,3,10, 4,9,9,9, 6,8,8,8, 11,10,10,7, 10,11,11,11]

orb = dest_to_orb(Dest)

tet0 = orb.Tetrahedra[0]
tet1 = orb.Tetrahedra[1]
tet2 = orb.Tetrahedra[2]

next_orb = CuspedOrbifold(two_to_three(orb.Tetrahedra,tet1,F3))

tet0 = next_orb.Tetrahedra[0]
tet1 = next_orb.Tetrahedra[1]

#show_triangulation(next_orb.Tetrahedra)


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

#show_triangulation(next_orb.Tetrahedra)

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

#print(next_orb.is_canonical)
# apparently this is canonical.
"""





"""
The following program takes the dest seqs from enum36success.txt and makes them into a list, Dest_Seqs,
where each element of the list is an integer list representing a dest seq. Then it turns them into
Cuspedorbifold objects, storing them in the list Orbs_up_to_36. We also sort those into two lists,
Canonical_Orbs and Non_Canonical_Orbs, according to whether or not the decomposition is canonical.
If all tilt sums are non-positive, we call it canonical. So, in Goerner et al's language, they're actually
only proto-canonical.
"""
"""
file1 = open("enum36success.txt","r")

L_out = file1.readlines()

int_strings = []
for i in range(100):
    int_strings.append(str(i))

Dest_Seqs = []

for word in L_out:
    if word[0] in int_strings and word[1] == ",":
        int_list = []
        for i in range(len(word)):
            if i == 0 or word[i-1] == "," or word[i-1] == " ":
                if i != len(word) - 1:
                    int_at_i = word[i]
                    if word[i+1] in int_strings:
                        int_at_i = int_at_i + word[i+1]
                    int_list.append(int(int_at_i))
        Dest_Seqs.append(int_list)

Orbs_up_to_36 = []
Canonical_Orbs = []
Non_Canonical_Orbs = [] 

for Dest in Dest_Seqs:
    orb = dest_to_orb(Dest)
    Orbs_up_to_36.append(orb)
    if orb.is_canonical is True:
        Canonical_Orbs.append(orb)
    else:
        Non_Canonical_Orbs.append(orb)

print(len(Orbs_up_to_36))
print(len(Canonical_Orbs))
print(len(Non_Canonical_Orbs))
for i in range(10):
    print(Non_Canonical_Orbs[i].DestSeq)
"""
