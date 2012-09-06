"""
A class to enumerate all elliptic curves over Q up to isomorphism,
given in short Weierstrass form, ordered by height.
"""

import time
import numpy as np
from copy import deepcopy
from sage.rings.integer_ring import IntegerRing, ZZ
from sage.functions.other import floor,ceil
from sage.combinat.cartesian_product import CartesianProduct
from sage.misc.misc import powerset, srange
from sage.schemes.elliptic_curves.constructor import EllipticCurve

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
            self._num_coeffs = ZZ(2)
            self._pows = (ZZ(3),ZZ(2))

    def __repr__(self):
        """
        Representation of self. Prints what family of curves is being considered,
         the model description and the height function on that family.
        
        EXAMPLES::

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: C
            Height iterator for elliptic curves over Q
            Family:             Short Weierstrass
            Model:              Y^2 = X^3 + A*X + B
            Height function:    H = min{|A|^3,|B|^2}
        """

        if self._family=="short_weierstrass":
            name = "Short Weierstrass"
            model = "Y^2 = X^3 + A*X + B"
            height_function = "min{|A|^3,|B|^2}"

        s = "Height iterator for elliptic curves over Q\n"
        s += "Family:             "+name+"\n"
        s += "Model:              "+model+"\n"
        s += "Height function:    H = "+height_function
        return s

    def _height_increment(self,coeffs):
        """
        Given a tuple of coefficients of a Weierstrass equation with
         a certain height, return the next largest permissable height
         and the range of coefficients that achieve that height.

        INPUT:
        
            - ``coeffs`` -- A list or tuple of coefficients of the same
              length as the number of coeffients in the model.

        OUTPUT:

            - A tuple of three entries consisting of:
              -- the smallest permissable height greater than the height
                 of the input coefficient list;
              -- a list of coeffients, some of which attain the above
                 height;
              -- a list of indices of which of the coefficients in the
                 above list achieve this height. The remaining entries
                 in the coefficient list indicate the maximum absolute
                 value that coefficient can attain without affecting the
                 curve's height.  

        EXAMPLES::

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: C._height_increment([2,2])
            (9, [2, 3], [1])
            sage: C._height_increment([3,7])
            (64, [4, 8], [0, 1]) 
        """
        I = range(len(coeffs))
        height_candidates = [(coeffs[i]+1)**(self._pows[i]) for i in I]
        next_height = min(height_candidates)
        index = []

        new_coeffs = list(coeffs)
        for i in I:
            if height_candidates[i]==next_height:
                new_coeffs[i] += 1
                index.append(i)
        return (next_height,new_coeffs,index)

    def next_height(self,N):
        """
        Return the next permissable height greater than or equal to N for
         curves in self's family.
        
        WARNING: This function my return a height for which only singular
                  curves exist. For example, in the short Weierstrass case
                  height 0 is permissable, as the curve Y^2 = X^3 (uniquely)
                  has height zero.

        INPUT:

            - ``N`` -- A non-negative integer

        OUTPUT:

            - A tuple consisting of three elements of the form
              (H, C, I) such that 
              H: The smallest height >= N
              C: A list of coefficients for curves of this height
              I: A list of indices indicating which of the above coefficients
              achieve this height. The remaining values in C  indicate the 
              max absolute value those coefficients are allowed to obtain
              without altering the height.

              For example, the tuple (4, [1, 2], [1]) for the short Weierstrass
              case denotes set of curves with height 4; these are all of the
              form Y^2 = X^3 + A*X + B, where B=2 and A ranges between -1 and 1.

        EXAMPLES::

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: C.next_height(4) 
            (4, [1, 2], [1])
            sage: C.next_height(60)
            (64, [4, 8], [0, 1])

            sage: C.heights(-100)
            Traceback (most recent call last):
            ...
            AssertionError: Input must be non-negative integer.
        """

        assert N>=0, "Input must be non_negative integer."

        coeffs = [ceil(N**(1/n))-1 for n in self._pows]
        height = max([coeffs[i]**(self._pows[i]) for i in range(self._num_coeffs)])
        return self._height_increment(coeffs)

    def heights(self,lowerbound,upperbound):
        """
        Return a list of permissable curve heights in the specified range
         (bounds inclusive), and for each height the equation coefficients
         that produce curves of that height.

        WARNING: This function my return heights for which only singular
                  curves exist. For example, in the short Weierstrass case
                  height 0 is permissable, as the curve Y^2 = X^3 (uniquely)
                  has height zero. 

        INPUT:

            - ``lowerbound`` -- Lower bound for the height range;
            - ``upperbound`` -- Upper bound for the height range. Heights
              returned are up to and including both bounds.

        OUTPUT:

            - A list of tuples, each consisting of three elements of the form
              (H, C, I) such that 
              H: The smallest height >= N
              C: A list of coefficients for curves of this height
              I: A list of indices indicating which of the above coefficients
              achieve this height. The remaining values in C  indicate the 
              max absolute value those coefficients are allowed to obtain
              without altering the height.

              For example, the tuple (4, [1, 2], [1]) for the short Weierstrass
              case denotes set of curves with height 4; these are all of the
              form Y^2 = X^3 + A*X + B, where B=2 and A ranges between -1 and 1.

        EXAMPLES::

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: C.heights(100,150) 
            [(121, [4, 11], [1]), (125, [5, 11], [0]), (144, [5, 12], [1])]

            sage: C.heights(-100,100)
            Traceback (most recent call last):
            ...
            AssertionError: Height upper bound must be greater than or equal to lower bound.
            sage: C.heights(150,100)
            Traceback (most recent call last):
            ...
            AssertionError: Height upper bound must be non-negative.
        """

        assert lowerbound>=0, "Height lower bound must be non-negative."
        assert upperbound>=lowerbound, "Height upper bound must be greater "\
                   +"than or equal to lower bound." 

        coeffs = [ceil(lowerbound**(1/n))-1 for n in self._pows]
        height = max([coeffs[i]**(self._pows[i]) for i in range(self._num_coeffs)])

        L = []
        while height <= upperbound:
            C = self._height_increment(coeffs)
            if C[0]>upperbound:
                break
            else:
                height = C[0]
                coeffs = C[1]
                L.append(C)
        return L

    def _coeffs_to_a_invariants(self,c):
        """
        Convert curve coefficients to a-invariants. This is family-specific.

        INPUT:

            - ``c`` -- The list of coefficients of the equation of a curve
              as per the family model description. See the __init__() method
              of this class for more info.

        OUTPUT:

            - A list of five integers corresponding to the a-invariants of
              the curve.

        EXAMPLES::
            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: C._coeffs_to_a_invariants([4,9])
            [0, 0, 0, 4, 9]
        """
        if self._family == "short_weierstrass":
            C = [0,0,0,c[0],c[1]]
        elif self._family == "rank_one":
            C = [0,c[0],c[1],c[2],0]
        elif self._family == "rank_two":
            a2,a4,b4,a6 = c[0],c[1],c[2],c[3]
            I = 3*(a4**2)+b4**2-3*a2*a6
            J = -27/4*(a2**2)*a4**2+18*(a4**2)*b4 - 2*(b4**3) \
                + 9*a2*b4*a6 - 27*(a6**2)
            if J.denominator()==4:
                I = I*16
                J = J*64
            C = [0,0,0,-27*I,-27*J]
        elif self._family == "two_torsion":
            C = [0,c[0],0,c[1],0]
        elif self._family == "three_torsion":
            C = [c[0],0,c[1],0,0]
        elif self._family == "full_weierstrass":
            raise NotImplementedError("Not yet implemented.")
        elif self._family == "F_1(2)":
            C = [0,c[0]**2-c[1]-c[2],0,c[1]*c[2],0]

        return C

    def _is_singular(self,C):
        """
        Tests if the a-invariants in 5-tuple C specify a singular
         elliptic curve.

        INPUT:

            - ``C`` -- A 5-tuple/list of a-invariants of a potential
              elliptic curve over Q

        OUTPUT:

            - True or False

        EXAMPLES::
            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: C.is_singular([0,0,0,3,2])
            False
            sage:EllipticCurve([0,0,0,3,2])
            Elliptic Curve defined by y^2 = x^3 + 3*x + 2 over Rational Field
            sage: C.is_singular([0,0,0,-3,2])
            True
            sage: EllipticCurve([0,0,0,-3,2])
            Traceback (most recent call last):
            ...
            ArithmeticError: Invariants [0, 0, 0, -3, 2] define a singular curve.
        """
        a1 = C[0]
        a2 = C[1]
        a3 = C[2]
        a4 = C[3]
        a6 = C[4]

        b2 = a1**2 + 4*a2
        b4 = 2*a4 + a1*a3
        b6 = a3**2 + 4*a6
        b8 = (a1**2)*a6 + 4*a2*a6 - a1*a3*a4 + a2*(a3**2) - a4**2

        Delta = -(b2**2)*b8 - 8*(b4**3) - 27*(b6**2) + 9*b2*b4*b6
        return Delta==0

    def _coefficients_from_height(self,height,coeffs,index):
        """
        Returns a list of tuples of a-invariants of all curves
         of the specified height.

        INPUT:

            - ``height`` -- A permissable curve height
            - ``coeffs`` -- A tuple of coeffients as per the family model
            - ``index``  -- The index of coefficients in the above tuple
              that achieve this height. See the documentation for 
              _height_increment() for a more thorough descripltion of this
              format.

        OUTPUT:

            - A list of 2-tuples, each consisting of the given height,
              followed by a tuple of a-invariants of a curve of that height.

        EXAMPLES:

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: B = C.heights(4,4)[0]; B
            (4, [1, 2], [1])
            sage: L = C.coefficients_from_height(*B); L
            (4, [0, 0, 0, -1, -2]), (4, [0, 0, 0, -1, 2]), (4, [0, 0, 0, 0, -2]), \
            (4, [0, 0, 0, 0, 2]), (4, [0, 0, 0, 1, -2]), (4, [0, 0, 0, 1, 2])]
        """
        # Produce list of all coefficient tuples with given height
        L = []
        for S in list(powerset(index))[1:]:
            B = []
            for j in range(len(coeffs)):
                if j in S:
                    B.append([-coeffs[j],coeffs[j]])
                elif j in index:
                    B.append(srange(-coeffs[j]+1,coeffs[j]))
                else:
                    B.append(srange(-coeffs[j],coeffs[j]+1))
            C = CartesianProduct(*B).list()
            for c in C:
                L.append(c)

        # Convert coefficient tuples to a-invariants
        L2 = []
        for c in L:
            C = self._coeffs_to_a_invariants(c)
            if not self._is_singular(C):
                L2.append((height,C))
        return L2

    def _coefficients_from_height_list(self,coefficient_list):
        """
        Return all height/a-invariant tuples of elliptic curves from a
         list of curve height/coefficient/index tuples.

        INPUT:

            - ``coefficient_list`` -- A list of height/coefficient/index
              tuples. See the documentation for _height_increment() for
              a description of this format.

        OUTPUT:

            - A list of tuples, each consisting of a height and a tuple of
              a-invariants defining an elliptic curve over Q of that height.
              The list will be ordered by increasing height.

        EXAMPLES:

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: B = C.heights(1,4)[0]; B
            [(1, [1, 1], [0, 1]), (4, [1, 2], [1])]
            sage: L = C.coefficients_from_height_list(B)
            sage: for ell in L: print(ell)
            ....: 
            (1, [0, 0, 0, -1, 0])
            (1, [0, 0, 0, 1, 0])
            (1, [0, 0, 0, 0, -1])
            (1, [0, 0, 0, 0, 1])
            (1, [0, 0, 0, -1, -1])
            (1, [0, 0, 0, -1, 1])
            (1, [0, 0, 0, 1, -1])
            (1, [0, 0, 0, 1, 1])
            (4, [0, 0, 0, -1, -2])
            (4, [0, 0, 0, -1, 2])
            (4, [0, 0, 0, 0, -2])
            (4, [0, 0, 0, 0, 2])
            (4, [0, 0, 0, 1, -2])
            (4, [0, 0, 0, 1, 2])
        """
        L2 = []
        for C in coefficient_list:
            L2 += self._coefficients_from_height(*C)
        return L2

    def coefficients_over_height_range(self,lowerbound,upperbound,\
                                       savefile=None,return_data=True):
        """
        Return all a-invariant tuples of elliptic curves over a given
         height range, bounds inclusive.

        INPUT:

            - ``lowerbound``  -- The lower height bound
            - ``upperbound``  -- The upper height bound
            - ``savefile``    -- String: the name of a file to which the
              output will be saved. Each line of the save file describes
              a single curve, and consists of six tab-separated integers;
              the first is the height of the curve, and the following five
              its a-invariants.
            - ``return_data`` -- (Default: True) If False, the computed
              data will not be returned.

        OUTPUT:

            - A list of tuples, each consisting of a height and a tuple of
              a-invariants defining an elliptic curve over Q of that height.
              The list will be ordered by increasing height.

        EXAMPLES:

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: L = C.coefficients_over_height_range(0,4)
            sage: for ell in L: print(ell)
            ....: 
            (1, [0, 0, 0, -1, 0])
            (1, [0, 0, 0, 1, 0])
            (1, [0, 0, 0, 0, -1])
            (1, [0, 0, 0, 0, 1])
            (1, [0, 0, 0, -1, -1])
            (1, [0, 0, 0, -1, 1])
            (1, [0, 0, 0, 1, -1])
            (1, [0, 0, 0, 1, 1])
            (4, [0, 0, 0, -1, -2])
            (4, [0, 0, 0, -1, 2])
            (4, [0, 0, 0, 0, -2])
            (4, [0, 0, 0, 0, 2])
            (4, [0, 0, 0, 1, -2])
            (4, [0, 0, 0, 1, 2])
        """
        H = self.heights(lowerbound,upperbound)
        L = self._coefficients_from_height_list(H)

        # Save data to file
        if savefile is not None:
            out_file  = open(savefile,"w")
            for C in L:
                out_file.write(str(C[0])+"\t")
                for a in C[1]:
                    out_file.write(str(a)+"\t")
                out_file.write("\n")

        if return_data:
            return L

    def _set_magma_class_group_bounds(proof = True):
        """
        Auxiliary function; set class group computation method
         in Magma (used in computing 3-Selmer group size)

        INPUT:

            - ``proof`` -- True or False; If False, the Generalized
              Riemann Hypothesis will not be assumed in the computing
              of class group bounds in Magma (and thus will be
              slower); if True, GRH will be assumed and computation
              will eb quicker.

        EXAMPLES:
        """
        if proof == False:
            magma.eval('SetClassGroupBounds("GRH")')
        else:
            magma.eval('SetClassGroupBounds("PARI")')


    def rank(self,curves,output_filename="rank_output.txt",\
             problems_filename="rank_problems.txt",\
             return_data=True,use_database=True,use_only_mwrank=False,\
             print_timing=True):
        r"""
        Compute the algebraic rank for a list of curves ordered by height.

        INPUT:
 
            - ``curves``            -- A list of height/a-invariant tuples of
              curves, as returned by the coefficients_over_height_range() method
            - ``output_filename``   -- String, the name of the file to which the
              output will be saved. Each line of the save file describes a
              single curve, and consists of seven tab-separated integers:
              the first is the height of the curve; the following five are
              the curve's a-invariants, and the final integer is the curve's rank.
            - ``problems_filename`` -- String: the file name to which problem
              curves will be written. These are curves for which rank could not
              be computed. Write format is the same as above, except rank is
              omitted at the end.
            - ``return_data``       -- (Default True): If set to False, the data
              is not returned at the end of computation; only written to file.
            - ``use_database``      -- (Default True): If set to False, the 
              Cremona database is not used to look up curve ranks.
            - ``use_only_mwrank``   -- (Default False): If set to True, will not
              try to use analytic methods to compute rank before consulting the
              Cremona database
            - ``print_timing``      -- (Default True): If set to False, wall time
              of total computation will not be printed.

        OUTPUT:

            - (only if return_data==True) A list consisting of two lists:
              The first is a list of triples, consisting of curve height, curve
              a-invariants, and computed rank. The second is a list of curves
              for which rank could not be (provably) computed; each entry of this
              list is just a pair consisting of height and a-invariants.

        EXAMPLES::

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: L = C.coefficients_over_height_range(0,4)
            sage: R = C.rank(L,return_data=true,print_timing=False)
            sage: R[1]
            []
            sage: for r in R[0]: print(r)
            ....: 
            (1, [0, 0, 0, -1, 0], 0)
            (1, [0, 0, 0, 1, 0], 0)
            (1, [0, 0, 0, 0, -1], 0)
            (1, [0, 0, 0, 0, 1], 0)
            (1, [0, 0, 0, -1, -1], 0)
            (1, [0, 0, 0, -1, 1], 1)
            (1, [0, 0, 0, 1, -1], 1)
            (1, [0, 0, 0, 1, 1], 1)
            (4, [0, 0, 0, -1, -2], 1)
            (4, [0, 0, 0, -1, 2], 0)
            (4, [0, 0, 0, 0, -2], 1)
            (4, [0, 0, 0, 0, 2], 1)
            (4, [0, 0, 0, 1, -2], 0)
            (4, [0, 0, 0, 1, 2], 0)
        """
        if print_timing:
            t = time.time()

        out_file  = open(output_filename,"w")
        prob_file = open(problems_filename,"w")
        if return_data:
            output = []
            problems = []
        for C in curves:
            # Attempt to compute rank and write curve+rank to file
            try:
                E = EllipticCurve(C[1])
                d = E.rank(use_database=use_database,only_use_mwrank=use_only_mwrank)
                out_file.write(str(C[0])+"\t")
                for a in C[1]:
                    out_file.write(str(a)+"\t")
                out_file.write(str(d)+"\n")
                out_file.flush()

                if return_data:
                    output.append((C[0],C[1],d))
            # Write to problem file if fail
            except:
                prob_file.write(str(C[0])+"\t")
                for a in C[1]:
                    prob_file.write(str(a)+"\t")
                prob_file.write("\n")
                prob_file.flush()

                if return_data:
                    problems.append(C)
        out_file.close()
        prob_file.close()

        if print_timing:
            print(time.time()-t)
        if return_data:
            return output,problems

    def averaged_data(self,input, outfile="average_data.txt", return_data=True):
        """
        INPUT:
            - ``input`` -- Either: A list where each element is
              (H, (a1,a2,a3,a4,a6), d), and
              H: curve height
              a1..a6: list of a-invariants of elliptic curve giving that height
              d: the type of data that's being averaged, e.g., rank, 2-Selmer,
              Or: String name of a file of array of data, where each line is
              H, a1, a2, a3, a4, a6, d
              where the elements are as above

            - ``outfile`` -- String: The name of the output file
            - ``return_data`` -- (Default: True) If set to False, the computed
              data is not returned.

        OUTPUT:
            - Writes the averaged data to file, where each line is of the form
              H, a
              where H is height, and a the average datum up to that height
            - (only if return_data==True) returns a list of tuples of the form
              (H, a)
              where H is height, and a the average datum up to that height

        EXAMPLES::

            sage: from ? import *
            sage: C = CurveEnumerator(model="short_weierstrass")
            sage: L = C.coefficients_over_height_range(0,4)
            sage: R = C.rank(L,return_data=true,print_timing=False)
            sage: A = C.averaged_data(R,return_data=True)
            sage: for c in A: print(c)
            ...:
        """
        # Load, format data
        if isinstance(input,str):
            data = np.loadtxt(input)
            X = data[:,0]
            V = data[:,-1]
        elif isinstance(input,list):
            X = np.array([C[0] for C in input])
            V = np.array([C[-1] for C in input])
        else:
            raise IOError("Input must either be a string or a list")

        #Compute running average
        N = np.arange(X.shape[0],dtype=np.float64) + 1
        Y = np.cumsum(V)/N

        # Retain only lines with new heights
        I = X[:-1] != X[1:]
        I = np.append(I,True)

        # Save, return output
        Z = np.vstack([X[I],Y[I,:].T]).T
        np.savetxt(outfile, Z)
        if return_data:
           return [(C[0],C[1]) for C in Z]


    def data_by_height2(L,inv="two_selmer",proof=True,return_data=False,\
                        output_filename="output.txt",problems_filename="problems.txt"):

        def _invariant(inv):
            if inv == "rank":
                def f(E): return E.rank(use_database=True, only_use_mwrank = False)
            elif inv == "sha_order":
                def f(E): return E.sha().an_numerical(prec=14)
            elif inv == "two_selmer_size":
                def f(E): return 2^(E.selmer_rank())
            elif inv == "two_selmer_rank":
                def f(E): return E.selmer_rank()
            elif inv == "reduced_two_selmer_size":
                def f(E): return 2^(E.selmer_rank()-E.two_torsion_rank())
            elif inv == "reduced_two_selmer_rank":
                def f(E): return E.selmer_rank()-E.two_torsion_rank()
            elif inv == "three_selmer_size":
                def f(E): return 3^(E.three_selmer_rank())
            elif inv == "three_selmer_rank":
                def f(E): return E.three_selmer_rank()
            elif inv == "reduced_three_selmer_size":
                def f(E): return 3^(E.three_selmer_rank()-valuation(E.torsion_order(),3))
            elif inv == "reduced_three_selmer_rank":
                def f(E): return E.three_selmer_rank()-valuation(E.torsion_order(),3)
            elif inv == "four_selmer_size":
                def f(E): return four_selmer_size(E)
            elif inv == "two_torsion_size":
                def f(E): return 2^(E.two_torsion_rank())
            elif inv == "two_torsion_rank":
                def f(E): return E.two_torsion_rank()
            elif inv == "three_torsion_size":
                def f(E): return 3^(valuation(E.torsion_order(),3))
            elif inv == "three_torsion_rank":
                def f(E): return valuation(E.torsion_order(),3)
            else: raise NotImplementedError('Invariant not yet Implemented')
            return f

        if inv in ["three_selmer_rank","three_selmer_size",\
                         "reduced_three_selmer_size",\
                         "reduced_three_selmer_rank"]:
            set_magma_class_group_bounds(proof)

        f =  _invariant(inv)

        t=time.time()
        out_file  = open(output_filename,"w")
        prob_file = open(problems_filename,"w")
        if return_data:
            output = []
            problems = []
        for C in L:
            try:
                E = EllipticCurve(C[1])
                d = f(E)
                out_file.write(str(C[0])+"\t")
                for a in C[1]:
                    out_file.write(str(a)+"\t")
                out_file.write(str(d)+"\n")
                if return_data:
                    output.append((C[0],C[1],f(E)))
            except:
                prob_file.write(str(C[0])+"\t")
                for a in C[1]:
                    prob_file.write(str(a)+"\t")
                prob_file.write("\n")
                if return_data:
                    problems.append(C)
            sys.stdout.flush()
        out_file.close()
        prob_file.close()
        print(time.time()-t)
        if return_data:
            return output,problems
