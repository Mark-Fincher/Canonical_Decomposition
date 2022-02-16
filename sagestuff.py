# sage stuff
#from sage.all import*
import sage.all
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.real_mpfi import RealIntervalField
from Exact_Arithmetic import*

RIF = RealIntervalField()
CIF = ComplexIntervalField()

def evaluate(self):
	"""
	Return an interval containing the true value.
	"""
	return sum([ RIF(r).sqrt() * RIF(c.numerator) / RIF(c.denominator) for r, c in self._entries])

def complex_evaluate(self):
	"""
	Returns a complex interval returning the true value.
	"""
	return CIF(self.real.evaluate(), self.imag.evaluate())

SquareRootCombination.evaluate = evaluate
ComplexSquareRootCombination.evaluate = complex_evaluate

"""
z = ComplexSquareRootCombination(SquareRootCombination.Zero(),SquareRootCombination.One())
print(z.imag.evaluate()) 
"""
