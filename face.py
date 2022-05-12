class Face:
	def __init__(self):
		self.Index = -1
		self.Corners = []

	def __repr__(self):
		if self.Index > -1:
			return ('f' + str(self.Index))        
		else:
			return '< floating face' + str(id(self)) + ' >'

	# Return a Corner whose two_subsimplex is not glued to None.
	def get_glued_corner(self):
		for corner in self.Corners:
			if corner.Tetrahedron.Neighbor[corner.Subsimplex] is not None:
				return corner