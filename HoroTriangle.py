# Horotriangle stuff from Goerner et al.

from fractions import Fraction
import sys
import math
from simplex import*
from Exact_Arithmetic import*

#---------T3M helper code---------

FacesAnticlockwiseAroundVertices = {
    V0 : (F1, F2, F3),
    V1 : (F0, F3, F2), 
    V2 : (F0, F1, F3),
    V3 : (F0, F2, F1)}

def glued_to(tetrahedron, face):
    """
    Returns (other tet, other face).
    """
    return tetrahedron.Neighbor[face], tetrahedron.Gluing[face].image(face)

def tets_and_vertices_of_cusp(cusp):
    return [(corner.Tetrahedron, corner.Subsimplex) for corner in cusp.Corners]


#--------HoroTriangle-----------

class HoroTriangle:
    """
    A horosphere cross section in the corner of an ideal tetrahedron.
    The sides of the triangle correspond to faces of the tetrahedron.
    """
    def __init__(self, tet, vertex, known_side, length_of_side):
        sides = FacesAnticlockwiseAroundVertices[vertex]
        left, center, right = HoroTriangle._make_middle(sides, known_side)
        z_l = tet.edge_params[left & center]
        z_r = tet.edge_params[center & right]
        L = length_of_side
        self.lengths = {center:L, left:abs(z_l)*L, right:L/abs(z_r)}
        a, b, c = self.lengths.values()

        real_arithmetic_type = type(z_l.real)
        Zero = real_arithmetic_type.Zero()
        One = real_arithmetic_type.One()
        Two = One + One
        Four = Two + Two
        self.area = L * L * z_l.imag / Two
        # Below is the usual formula for circumradius combined with 
        # Heron's formula.
        
        # 1/6/2022 Mark. Flat tets will have vertex horotriangles with area 0,
        # so undefined circumradius. In that case, set circumradius to None.
        # We don't care about making sense of tilt for flat tets. Canonize will
        # try to either cancel flat tets, or just ignore them until they can be
        # cancelled.
        if self.area != Zero:
            self.circumradius = (a * b * c /
                                 (Four * self.area))
        else:
            self.circumradius = None

    def rescale(self, t):
        "Rescales the triangle by a Euclidean dilation"
        for face, length in self.lengths.items():
            self.lengths[face] = t*length
        if self.circumradius is not None:
            self.circumradius = t*self.circumradius
        self.area = t*t*self.area

    @staticmethod
    # Changed by Mark, 7/12/2021, because tuple parameters aren't supported in python 3
    def _make_middle(triple, x):
        "Cyclically rotate (a,b,c) so that x is the middle entry"
        a,b,c = triple
        if x == a:
            return (c,a,b)
        elif x == b:
            return (a,b,c)
        elif x == c:
            return (b,c,a)



