"""
Cusped orbifold class. Like the t3m triangulation class, except we only include stuff we need,
and our triangulations are of orbifolds. Much of this is copied from Goerner et al, and that is
built on snappy.
"""
class CuspedOrbifold:
	
	def __init__(self,tets_list):
		self.Tetrahedra = tets_list
		for T in self.Tetrahedra:
            T.horotriangles = {V0:None, V1:None, V2:None, V3:None}
        self.add_cusp_cross_sections()

    