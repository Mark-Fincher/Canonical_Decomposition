from The_Algorithm import*
import json
import copy

"""
The following program takes the dest seqs from enum36success.txt and makes them into a list, Dest_Seqs,
where each element of the list is an integer list representing a dest seq. Then it turns them into
CuspedOrbifold objects, storing them in the list Orbs_up_to_36. We also sort those into two lists,
Canonical_Orbs and Non_Canonical_Orbs, according to whether or not the decomposition is canonical.
If all tilt sums are non-positive, we call it canonical. So, in Goerner et al's language, they're actually
only proto-canonical.

At this point, each orbifold is written in terms of the triangulation by regular ideal tetrahedra, which is
the first triangulation we get from the dest seq. The Non_Canonical_Orbs list consists of those for which this
is not the canonical decomposition. The last part of this program takes those orbifolds and attempts to find their
canonical decompositions using the "canonize" function in The_Algorithm.py. So far, this function only uses 2-3 moves.
If "canonize" succeeds in getting the canonical decomp, we put the new CuspedOrbifold in the list newly_canonical.
If it fails, then the orbifold is unchanged, and we put it in the list stuck_orbs.

Some simple output statistics are:

There are 175 total orbifolds
Canonical_Orbs has 124
Non_Canonical_Orbs has the remaining 51
Of those 51, 24 can be made canonical by canonize and the remaining 27 get stuck 

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

"""
OrbDictionary = {(1,0):Dest_Seqs[0]}

max_tets = 1

for dest in Dest_Seqs:
    if max_tets < len(dest)//4:
        max_tets = len(dest)//4

for n in range(max_tets):
    n = n + 1
    n_tet_Dest_seqs = []
    for dest in Dest_Seqs:
        if len(dest)//4 == n:
            n_tet_Dest_seqs.append(dest)
    sorted_n_tet_Dest_seqs = sorted(n_tet_Dest_seqs)
    for k in range(len(sorted_n_tet_Dest_seqs)):
        OrbDictionary[(n,k)] = sorted_n_tet_Dest_seqs[k]

OrbDictCopy = copy.deepcopy(OrbDictionary)

with open("OrbDictionary.json", "w") as write_file:
    keyz = OrbDictionary.keys()
    valz = OrbDictionary.values()
    json.dump({str(k):OrbDictionary[k] for k in keyz}, write_file)


with open("OrbDictionary.json", "r") as read_file:
    OrbDictionary = json.load(read_file)
    keyz = OrbDictionary.keys()
    OrbDictionary = {eval(k):OrbDictionary[k] for k in keyz}

print(OrbDictionary)
"""


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
print(len(Canonical_Orbs))
print(len(Non_Canonical_Orbs))

newly_canonical = []
stuck_orbs = []
for orb in Non_Canonical_Orbs:
    canonical_orb = canonize(orb)
    if canonical_orb is None:
        stuck_orbs.append(orb)
    else:
        newly_canonical.append(canonical_orb)

print(len(newly_canonical))
print(len(stuck_orbs))

"""
already_canonical_dest_seqs = []
for orb in Canonical_Orbs:
    already_canonical_dest_seqs.append(orb.DestSeq)

with open("already_canonical_dest_seqs.json", "w") as write_file:
    json.dump(already_canonical_dest_seqs, write_file)

newly_canonical_dest_seqs = []
for orb in newly_canonical:
    newly_canonical_dest_seqs.append(orb.DestSeq)

with open("newly_canonical_dest_seqs.json", "w") as write_file:
    json.dump(newly_canonical_dest_seqs, write_file)

stuck_dest_seqs = []
for orb in stuck_orbs:
    stuck_dest_seqs.append(orb.DestSeq)

with open("stuck_dest_seqs.json", "w") as write_file:
    json.dump(stuck_dest_seqs, write_file)
"""

