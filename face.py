class Face:
	def __init__(self):
		self.Index = -1
		self.Corners = []

	def __repr__(self):
		if self.Index > -1:
			return ('f' + str(self.Index))        
		else:
			return '< floating face' + str(id(self)) + ' >'