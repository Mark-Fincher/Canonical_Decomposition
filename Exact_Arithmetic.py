"""
This file is copied from part of Goerner et al's canonical_o3.py. I don't need the interval arithmetic,
so I comment out the parts relevant to that, including importing sage. This part of the code doesn't
use snappy, so I also comment that out.

This program creates two classes, SquareRootCombination and ComplexSquareRootCombination. Each is supposed to represent
a number which possibly has square roots, so wouldn't be stored in an exact form normally. In each class the definitions
of addition, multiplication, etc. are changed to be what they should be for that class. So, for example, if a and b 
are both objects in the ComplexSquareRootCombination class, running a*b will not produce an error (you would expect python
to not understand what it means to multiply things which aren't numbers), but will instead give the correct 
ComplexSquareRootCombination object. 
"""




# import snappy
# import snappy.snap.t3mlite as t3m
# from snappy.snap.t3mlite.simplex import *
import copy

#from sage.rings.complex_interval_field import ComplexIntervalField
#from sage.rings.real_mpfi import RealIntervalField

from fractions import Fraction
import sys
import math

#RIF = RealIntervalField()
#CIF = ComplexIntervalField()

class SquareRootCombination:
    """
    Represents a Q-linear combination of sqrt's of distinct positive
    square-free integers.
    Internally, c_1 * sqrt(r_1) + ... + c_n * sqrt(r_n) is represented by an
    array [(r_1, c_1), ..., (r_n, c_n)] such that r_i are ascending. The r_i
    are integers and the c_i integers or Fraction's.
    """

    @staticmethod
    def square_free(x):
        """
        Returns pair (t, y) such that x = t * y^2 and t is square-free.
        """

        y = 1
        i = 2
        # Not an efficient algorithm, but good enough for our purposes.
        while i * i <= x:
            
            if x % (i * i) == 0:
                x /= i * i
                y *= i
            else:
                i += 1
        
        return x, y

    def __init__(self, entries = []):
        """
        Constructs a Q-linear combination of square roots by
        normalizing the given entries [(root, coefficient), ...]
        """
        
        d = {}
        for root, coefficient in entries:
            square_free, extra_coefficient = (
                SquareRootCombination.square_free(root))
            d[square_free] = (
                d.get(square_free, 0) + coefficient * extra_coefficient)
        
        self._entries = sorted([(k,v) for (k,v) in d.items() if v])

    @staticmethod
    def Zero():
        return SquareRootCombination()

    @staticmethod
    def One():
        return SquareRootCombination([(1,1)])

    @staticmethod
    def Two():
        return SquareRootCombination([(1,2)])

    @staticmethod
    def Four():
        return SquareRootCombination([(1,4)])

    @staticmethod
    def SqrtThree():
        return SquareRootCombination([(3,1)])

    @staticmethod
    def guess_from_float(root, f):
        """
        Given an integer root and an object that can be converted to float f,
        tries to guess a representation of f as p/q * sqrt(root).
        """

        coeff = Fraction(float(f) / math.sqrt(root)).limit_denominator(10000)
        return SquareRootCombination([(root, coeff)])

    def __add__(self, other):
        return SquareRootCombination(self._entries + other._entries)

    def __neg__(self):
        return SquareRootCombination([(r, -c) for r, c in self._entries])

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return SquareRootCombination(
            [(r1*r2, c1*c2) for r1, c1 in self._entries
                            for r2, c2 in other._entries])

    """
    Chanded __div__ to __truediv__ because that's what python wants - Mark, 7/12/2021
    """
    def __truediv__(self, other):

        assert len(other._entries) > 0, "Division by zero not allowed"

        if len(other._entries) == 1:
            root, coefficient = other._entries[0]
            inv = SquareRootCombination([ (root,
                                           Fraction(1) / coefficient / root) ])
            return self * inv

        p = 2
        while True:
            divs = [(r, 2 * c) for r, c in other._entries if r % p == 0]
            if divs:
                f = other - SquareRootCombination(divs)
                numerator = self * f
                denominator = other * f
                return numerator / denominator
            p += 1
    
    def __str__(self):
        if not self._entries:
            return '0'
        return '+'.join(['(%s * sqrt(%d))' % (c, r) for (r, c) in self._entries])

    def __repr__(self):
        return 'SquareRootCombination(%r)' % self._entries

    def sqrt(self):
        err_msg = "Only square roots of rational numbers are supported"

        if len(self._entries) == 0:
            return SquareRootCombination([])
        assert len(self._entries) == 1, err_msg
        root, coefficient = self._entries[0]
        assert root == 1, err_msg
        assert coefficient > 0, err_msg
        if isinstance(coefficient, int):
            return SquareRootCombination([(coefficient, 1)])
        return (
            SquareRootCombination([(coefficient.numerator, 1)]) /
            SquareRootCombination([(coefficient.denominator, 1)]))

    def __eq__(self, other):
        return self._entries == other._entries

    #def evaluate(self):
        """
        Return an interval containing the true value.
        """
        #return sum([ RIF(r).sqrt() * RIF(c.numerator) / RIF(c.denominator)
         #            for r, c in self._entries])

class ComplexSquareRootCombination:
    """
    Represents a + b * i where a and b are Q-linear combinations of
    square roots of distinct square-free integers.

    This is implemented using SquareRootCombination objects stored
    under "real" and "imag".
    """


    def __init__(self, real, imag):
        """
        Constructs a + b * i given two SquareRootCombination objects.
        """

        self.real = real
        self.imag = imag

    @staticmethod
    def One():
        return ComplexSquareRootCombination(
            SquareRootCombination.One(), SquareRootCombination.Zero())

    def __repr__(self):
        return 'ComplexSquareRootCombination(%r, %r)' % (
            self.real, self.imag)

    def __abs__(self):
        return (self.real * self.real + self.imag * self.imag).sqrt()

    def __add__(self, other):
        return ComplexSquareRootCombination(
            self.real + other.real,
            self.imag + other.imag)

    def __neg__(self):
        return ComplexSquareRootCombination(
            -self.real, -self.imag)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return ComplexSquareRootCombination(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real)

    def conjugate(self):
        return ComplexSquareRootCombination(
            self.real, -self.imag)

    def __div__(self, other):
        otherConj = other.conjugate()
        denom = (other * otherConj).real
        num = self * otherConj
        return ComplexSquareRootCombination(
            num.real / denom, num.imag /denom)

    def __eq__(self, other):
        return self.real == other.real and self.imag == other.imag

    #def evaluate(self):
        """
        Returns a complex interval returning the true value.
        """

        #return CIF(self.real.evaluate(), self.imag.evaluate())
