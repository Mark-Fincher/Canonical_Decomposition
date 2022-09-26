# Canonical_Decomposition

This package is used to find canonical decompositions of hyperbolic 3-orbifolds in the commensurability class of the figure
eight knot complement. It will at some point be able to handle other non-compact finite-volume hyperbolic 3-orbifolds; the only real limitation currently 
is that the exact arithmetic program is tied to orbifolds in that commensurability class.

The main class is CuspedOrbifold. It represents a triangulation of a cusped hyperbolic 3-orbifold by hyperolic ideal tetrahedra. The tetrahedra are labelled with symmetry groups, the edges with integer labels, and a face of a tetrahedron is allowed to be glued to itself, or to another face. These things represent the orbifold structure of the orbifold which is represented by the triangulation. For more details, see a forthcoming preprint on my [website](https://sites.google.com/view/markfincher/home) soon.

# Credits

I started with SnapPy's [t3mlite](https://github.com/3-manifolds/SnapPy/tree/master/python/snap/t3mlite), then modified certain parts of it and added my 
own things. SnapPy is maintained by Marc Culler, Nathan Dunfield, and Matthias Goerner, and the original version was created by Jeff Weeks. I'm also 
using some code from the paper *A census of tetrahedral hyperbolic manifolds* by Fominykh, Garoufalidis, Goerner, Tarkaev, and Vesnin. Lastly, I use 
Sage for interval arithemtic.

# Example

The package must be used with sage. The file OrbDictionary.json contains the data of many finite index subgroups of PGL(2,O_3). Loading the file
gives a dictionary "OrbDictionary". OrbDictionary[(n,k)] is the "destination sequence" encoding a subgroup of index n in PGL(2,O_3) at lexicographic
position k among other destination sequences of the same index. For example, OrbDictionary[(1,0)] is PGL(2,O_3). In the below example, we turn it
into a CuspedOrbifold object, get its info, then check if it's proto-canonical or not (it is). 

```
from canonize import*
import json

with open("OrbDictionary.json", "r") as read_file:
    OrbDictionary = json.load(read_file)
    keys = OrbDictionary.keys()
    OrbDictionary = {eval(k):OrbDictionary[k] for k in keys}
    
orb = dest_to_orb(OrbDictionary[(1,0)])
orb.info()
print(orb.is_proto_canonical())
    
```

