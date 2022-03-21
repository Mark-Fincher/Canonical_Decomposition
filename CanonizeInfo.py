"""
canonize_info, for use in canonize_part2.py. Meant to be an attribute of a Tetrahedron.
"""
from simplex import*
from tetrahedron import*

class CanonizeInfo:
	def __init__(self):
		# Boolean
		self.part_of_coned_cell = None
		self.is_flat = None
		# face_status[face] == 0 if the face is opaque, 1 if the face is transparent,
		# 2 if the face is inside the coned polyhedron, or 3 if the face is glued to None.
		self.face_status = {F0:None,F1:None,F2:None,F3:None}

