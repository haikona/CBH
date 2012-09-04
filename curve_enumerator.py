"""
A class to enumerate all elliptic curves over Q up to isomorphism,
given in short Weierstrass form, ordered by height.
"""

from copy import deepcopy
from sage.rings.integer_ring import IntegerRing, ZZ

class CurveEnumerator():
    r"""
    The CurveEnumerator class will enumerate all elliptic curves
    over Q in a specified Weierstrass form up to isomorphism within
    a given height range, where height is a function of the Weierstrass
    equation of the curve.
    """

    def __init__(self, family="short_weierstrass"):
        """
        INPUT:


        EXAMPLES::

            sage:
        """
        
        self._family = family
        if family=="short_weierstrass":
            self._num_constants = ZZ(2)
            self._pows = ((ZZ(3),ZZ(2))

    def _height_increment(self,coeffs):
        """
        Given a tuple of coefficients of a Weierstrass equation with
         a certain height, return the next largest permissable height
         and the range of coefficients that achieve that height.

        INPUT:

        OUTPUT:

        EXAMPLES:: 
        """
        I = range(len(coeffs))
        height_candidates = [(coeffs[i]+1)**(self._pows[i]) for i in I]
        next_height = min(J)
        index = []

        new_coeffs = deepcopy(coeffs)
        for i in I:
            if height_candidates[i]==next_height:
                coeffs2[i] += 1
                index.append(i)
        return (next_height,new_coeffs,index)

    def heights(self,lowerbound=0,upperbound):
        """
        Return a list of permissable curve heights in the specified range,
         and for each height the equation coefficients that produce curves
         of that height.

        INPUT:

        OUTPUT:

        EXAMPLES:: 
        """
        
