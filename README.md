# Canonical_Decomposition

This package, still under construction, is used to find canonical decompositions of hyperbolic 3-orbifolds in the commensurability class of the figure
eight knot complement. It will at some point be able to handle other non-compact finite-volume hyperbolic 3-orbifolds; the only real limitation currently 
is that the exact arithemtic program is tied to orbifolds in that commensurability class.

The main class is CuspedOrbifold. It represents a fully ideal triangulation. The main difference between an ideal triangulation of an orbifold and
an ideal triangulation of a manifold is that tetrehdra in the orbifold triangulation may have some symmetries acting on them. These orbifold triangulations
are not necessarily triangulations of the underlying space, which is more standard.

Preprints concerning this project should appear on my [website](https://sites.google.com/view/markfincher/home) soon.

# Credits

I started with SnapPy's [t3mlite](https://github.com/3-manifolds/SnapPy/tree/master/python/snap/t3mlite), then modified certain parts of it and added my 
own things. SnapPy is maintained by Marc Culler, Nathan Dunfield, and Matthias Goerner, and the original version was created by Jeff Weeks. I'm also 
using some code from the paper *A census of tetrahedral hyperbolic manifolds* by Fominykh, Garoufalidis, Goerner, Tarkaev, and Vesnin. Lastly, I use 
Sage for interval arithemtic.

# Example

The package must be used with sage. The file OrbDictionary.json contains the data of many finite index subgroups of PGL(2,O_3). Loading the file
gives a dictionary "OrbDictionary". OrbDictionary[(n,k)] is the "destination sequence" encoding a subgroup of index n in PGL(2,O_3) at lexicographic
position k among other destination sequences of the same index. For example, OrbDictionary[(1,0)] is PGL(2,O_3). In the below example, we turn it
into a CuspedOrbifold object, get its info, then check if it's already canonical or not. It is already canonical. 

```
from canonize import*
import json

with open("OrbDictionary.json", "r") as read_file:
    OrbDictionary = json.load(read_file)
    keyz = OrbDictionary.keys()
    OrbDictionary = {eval(k):OrbDictionary[k] for k in keyz}
    
orb = dest_to_orb(OrbDictionary[(1,0)])
orb.info()
print(orb.is_canonical)
    
```

