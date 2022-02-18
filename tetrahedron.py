#$Id: tetrahedron.py,v 1.2 2002/09/20 03:52:16 culler Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

# CREDITS. Most of this file is taken from tetrahedron.py of t3mlite (including the above).
# But I've added some orbifold things. For my own contributions, I comment
# ORBIFOLDS before it. I've also added some things from the ancillary files
# of the paper "A census of tetrahedral hyperbolic manifolds" by Evgeny Fominykh, 
# Stavros Garoufalidis, Matthias Goerner, Vladimir Tarkaev, and Andrei Vesnin. I comment
# FGGTV before those things.
# - Mark F. 

from simplex import *
import sys
from Exact_Arithmetic import*

class Tetrahedron:
    def __init__(self, name = ''):
        self.Index = -1
        self.Name = name
        self.Neighbor = {F0:None,F1:None,F2:None,F3:None}  # Tetrahedra
        self.Gluing   = {F0:None,F1:None,F2:None,F3:None}  # Permutations
        self.Class    = [None]*16             # list of equivalence classes
        self.Checked  = 0                     # flag
        
        # ORBIFOLDS
        self.Symmetries = []
        self.edge_group_labels = {E01:None,E23:None,E02:None,E13:None,E03:None,E12:None}

        # FGGTV
        self.edge_params = {E01:None,E23:None,E02:None,E13:None,E03:None,E12:None}
        self.horotriangles = {V0:None, V1:None, V2:None, V3:None}

    # FGGTV
    def tilt(self, v):
        "The tilt of the face of the tetrahedron opposite the vertex v."
        ans = SquareRootCombination.Zero()
        for w in ZeroSubsimplices:
            if v == w:
                c_w = SquareRootCombination.One()
            else:
                z = self.edge_params[v | w]
                c_w = -z.real/abs(z)
            R_w = self.horotriangles[w].circumradius

            ans += c_w*R_w
        return ans

    def __repr__(self):
        if self.Index != -1: 
            return ( 'tet'+ str(self.Index) )
        else:
            return '< floating tetrahedron ' + ' at ' + str(id(self)) + '>'

    def attach(self, two_subsimplex, tet, perm_data):
        if tet == None:
            self.Neighbor[two_subsimplex] = None
            self.Gluing[two_subsimplex] = None
        else:
            perm = Perm4(perm_data)
            self.Neighbor[two_subsimplex] = tet
            self.Gluing[two_subsimplex] = perm
            tet.Neighbor[perm.image(two_subsimplex)] = self
            tet.Gluing[perm.image(two_subsimplex)] = inv(self.Gluing[two_subsimplex])
# Reverse the orientation.  Vertices are relabelled by a transposition
# and gluings are adjusted.
#
    def reverse(self):
        transpo = Perm4((1,0,2,3))
        nhbr = self.Neighbor.copy()
        gluing = self.Gluing.copy()
        for two_subsimplex in TwoSubsimplices:
            relabeled = transpo.image(two_subsimplex)
            if not nhbr[two_subsimplex] == None:
                perm = (gluing[two_subsimplex]*transpo).tuple()
            else:
                perm = None
            self.attach(relabeled, nhbr[two_subsimplex], perm)

# Unglues and removes references to self from neighbor.
#
    def detach(self, two_subsimplex):
        neighbor = self.Neighbor[two_subsimplex]
        if neighbor == None:
            return
        neighbors_subsimplex = self.Gluing[two_subsimplex].image(two_subsimplex)
        self.Neighbor[two_subsimplex] = None
        self.Gluing[two_subsimplex] = None
        if (neighbor.Neighbor and 
                neighbor.Neighbor[neighbors_subsimplex] == self):
            neighbor.Neighbor[neighbors_subsimplex] = None
            neighbor.Gluing[neighbors_subsimplex] = None

    def erase(self):
        for two_subsimplex in TwoSubsimplices:
            self.detach(two_subsimplex)
        self.Index = -1
        self.Neighbor = None
        self.Gluing = None
        self.clear_Class()

    def clear_Class(self):
        self.Class    = [None]*16              # list of equivalence classes

    def info(self, out = sys.stdout):
        if len(self.Name) == 0:
            out.write(repr(self) + "\t%s\n" %
                      ([self.Neighbor.get(s) for s in TwoSubsimplices]))
        else:
            out.write(repr(self) + " ( " + self.Name + " )\n")
            out.write("\t%s\n" % ([self.Neighbor.get(s) for s in TwoSubsimplices]))

        out.write("\t%s\n" % ([self.Gluing.get(s) for s in TwoSubsimplices]))

        out.write("\tVertices: " + repr(self.Class[V0]) 
                                 + repr(self.Class[V1])
                                 + repr(self.Class[V2])
                                 + repr(self.Class[V3]) + '\n')

        if self.Index > -1:
            s = ""
            for edge in OneSubsimplices[:3]:
                s = (s + "%s : %-10s   " %
                     (SubsimplexName[edge], self.Class[edge]))
            out.write("\tEdges: " + s + '\n')
            s = ""
            for edge in OneSubsimplices[3:]:
                s = (s + "%s : %-10s   " % 
                     (SubsimplexName[edge], self.Class[edge]))
            out.write("\t       " + s + '\n')

    # Below added 7/12/99 by NMD

    def get_orientation_of_edge(self, a, b):
        return self.Class[a | b].orientation_with_respect_to(self, a, b)

    def fill_edge_params(self,z):
        One = ComplexSquareRootCombination.One()
        zp  = One / (One - z)
        zpp = (z - One) / z
        self.edge_params = {E01:z, E23:z, E02:zp, E13:zp, E03:zpp, E12:zpp}

    """
    ORBIFOLDS.
    The folllowing returns true if the face is glued to itself, false otherwise.
    """
    def face_glued_to_self(self,two_subsimplex):
        if self.Neighbor[two_subsimplex] == self and self.Gluing[two_subsimplex].image(two_subsimplex) == two_subsimplex:
            return True
        return False

    """
    ORBIFOLDS.
    My convention is that if the equivalence relation induced by the face identifications and
    symmetries of a tetrahedron imply that a certain face is glued to itself, then that should
    explicitly be the gluing data of that face. Note that a face could be glued to None but
    still be glued to itself by the induced equivalence relation. The following function
    makes sure this convention is followed for self.
    """
    def fix_glued_to_self(self):
        if len(self.Symmetries) == 1:
            return
        for face in TwoSubsimplices:
            if self.Neighbor[face] is not None and self.Neighbor[face] is not self:
                continue
            if self.face_glued_to_self(face):
                continue
            face_done = False
            for sym1 in self.Symmetries:
                if self.Neighbor[sym1.image(face)] is self:
                    for sym2 in self.Symmetries:
                        perm = sym2*self.Gluing[sym1.image(face)]*sym1 
                        if perm.image(face) == face:
                            self.detach(face)
                            self.attach(face,self,perm.tuple())
                            face_done = True
                            break
                if face_done:
                    break


    """
    ORBIFOLDS.
    The following returns true if there's a symmetry of the tet which rotates the face, false otherwise.
    """
    def face_rotation(self,two_subsimplex):
        for sym in self.Symmetries:
            if sym.image(two_subsimplex) == two_subsimplex and sym.tuple() != (0,1,2,3):
                return True
        return False

    """
    ORBIFOLDS.
    """
    def is_flat(self):
        z = self.edge_params[E01]
        if z.imag == SquareRootCombination.Zero():
            return True
        else:
            return False

    """
    ORBIFOLDS
    """
    def is_regular(self):
        a = SquareRootCombination([(1,Fraction(1/2))])
        b = SquareRootCombination([(Fraction(3,1),Fraction(1/2))])
        z = ComplexSquareRootCombination(a,b)
        if self.edge_params[E01] == z:
            return True
        else:
            return False

    """
    ORBIFOLDS
    """
    def true_glued(self,two_subsimplex):
        #return (neighbor,gluing map), and if two_subsimplex is glued to None
        #but there's a symmetry taking it to another face which is not glued to
        #None, use that data.
        if self.Neighbor[two_subsimplex] is not None:
            return (self.Neighbor[two_subsimplex],self.Gluing[two_subsimplex])
        for sym in self.Symmetries:
            if self.Neighbor[sym.image(two_subsimplex)] is not None:
                Neighbor = self.Neighbor[sym.image(two_subsimplex)]
                perm = self.Gluing[sym.image(two_subsimplex)]
                if self.face_glued_to_self(sym.image(two_subsimplex)):
                    return(Neighbor,sym.inverse()*perm*sym)
                else:
                    return(Neighbor,perm*sym)
        return (None,None)
