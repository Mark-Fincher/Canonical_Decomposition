# Horotriangle stuff from Goerner et al.

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



