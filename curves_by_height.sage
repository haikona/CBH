
# A more systematic attempt at code for generating statistics of invariants
# of elliptic curves ordered by height.

from sage.interfaces.all import magma
import time
import numpy as np


def height_iterator(M,N,model):
    if model == "short_weierstrass":
       	return height_iterator_short_weierstrass(M,N)
    elif model == "rank_one":
        return height_iterator_rank_one(M,N)
    elif model == "rank_two":
        return height_iterator_rank_two(M,N)
    elif model == "two_torsion":
        return height_iterator_two_torsion(M,N)
    elif model == "three_torsion":
        return height_iterator_three_torsion(M,N)
    elif model == "full_weierstrass":
        return height_iterator_full_weierstrass(M,N)
    else: raise IOError("Please enter recognized Weierstrass family of curves.")

def height_iterator_short_weierstrass(M,N):
    L = []

    a = floor(M^(1/3))
    b = floor(M^(1/2))
    H = max([a^3,b^2])
    while H <= N:
        a2 = (a+1)^3
        b2 = (b+1)^2
        if min([a2,b2]) > N:
            break
        if a2 < b2:
            a += 1
            H = a2
            bbound = floor(H^(1/2))
            L.append((H,[a,bbound],[0]))
        elif b2 < a2:
            b += 1
            H = b2
            abound = floor(H^(1/3))
            L.append((H,[abound,b],[1]))
        elif a2 == b2:
            a += 1
            b += 1
            H = a2
            L.append((H,[a,b],[0,1]))
    return L

def height_iterator_rank_one(M,N):
    L = []

    x2 = floor(M^(1/6))
    x3 = floor(M^(1/4))
    x4 = floor((M)^(1/3))

    H = max([x2^6,x3^4,x4^3])
    while H <= N:
        y2 = (x2+1)^6
        y3 = (x3+1)^4
        y4 = (x4+1)^3
        yl = [y2,y3,y4]
        if min([y2,y3,y4]) > N:
            break
        if y2 < min(y3,y4):
            x2 += 1
            H = y2
            x3bound = floor((H)^(1/4))
            x4bound = floor((H)^(1/3))
            L.append((H,[x2,x3bound,x4bound],[0]))

        if y3 < min(y2,y4):
            x3 += 1
            H = y3
            x2bound = floor((H)^(1/6))
            x4bound = floor((H)^(1/3))
            L.append((H,[x2bound,x3,x4bound],[1]))

        if y4 < min(y2,y3):
            x4 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            x3bound = floor((H)^(1/4))
            L.append((H,[x2bound,x3bound,x4],[2]))
            
        if y2 == y3 and y2 < y4:
            x2 += 1
            x3 += 1
            H = y2
            x4bound = floor((H)^(1/3))
            L.append((H,[x2,x3,x4bound],[0,1]))
            
        if y2 == y4 and y2 < y3:
            x2 += 1
            x4 += 1
            H = y2
            x3bound = floor((H)^(1/4))
            L.append((H,[x2,x3bound,x4],[0,2]))
           
        if y3 == y4 and y3 < y2:
            x3 += 1
            x4 += 1
            H = y2
            x2bound = floor((H)^(1/6))
            L.append((H,[x2bound,x2,x4],[1,2]))

        if y2 == y3 and y2 == y4:
            x2 += 1
            x3 += 1
            x4 += 1
            H = y2
            L.append((H,[x2,x3,x4],[0,1,2]))
    return L

def height_iterator_rank_two(M,N):
    L = []
    x2 = floor(M^(1/6))
    x4 = floor((M)^(1/3))
    xp4 = floor((M)^(1/3))
    x6 = floor((M)^(1/2))
    H = max([x2^6,x4^3,xp4^3,x6^2])
    while H <= N:
        y2 = (x2+1)^6        
        y4 = (x4+1)^3
        yp4 = (xp4+1)^3
        y6 = (x6+1)^2        
        yl = [y2,y4,yp4,y6]
        if min(yl) > N:
            break

        if y2 < min(y4,yp4,y6):
            x2 += 1
            H = y2
            x4bound = floor((H)^(1/3))
            xp4bound = floor((H)^(1/3))
            x6bound = floor((H)^(1/2))
            L.append((H,(x2,x4bound,xp4bound,x6bound),[0]))

        if y4 < min(y2,yp4,y6):
            x4 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            xp4bound = floor((H)^(1/3))
            x6bound = floor((H)^(1/2))
            L.append((H,(x2bound,x4,xp4bound,x6bound),[1]))

        if yp4 < min(y2,y4,y6):
            xp4 += 1
            H = yp4
            x2bound = floor((H)^(1/6))
            x4bound = floor((H)^(1/3))
            x6bound = floor((H)^(1/2))
            L.append((H,(x2bound,x4bound,xp4,x6bound),[2]))

        if y6 < min(y2,y4,yp4):
            x6 += 1
            H = y6
            x2bound = floor((H)^(1/6))
            x4bound = floor((H)^(1/3))
            xp4bound = floor((H)^(1/3))
            L.append((H,(x2bound,x4bound,xp4bound,x6),[3]))

        if y2 == y4 and y2 < min(yp4,y6):
            x2 += 1
            x4 += 1
            H = y2
            xp4bound = floor((H)^(1/3))
            x6bound = floor((H)^(1/2))
            L.append((H,(x2,x4,xp4bound,x6bound),[0,1]))

        if y2 == yp4 and y2 < min(y4,y6):
            x2 += 1
            xp4 += 1
            H = y2
            x4bound = floor((H)^(1/3))
            x6bound = floor((H)^(1/2))
            L.append((H,(x2,x4bound,xp4,x6bound),[0,2]))

        if y2 == y6 and y2 < min(y4,yp4):
            x2 += 1
            x6 += 1
            H = y2
            x4bound = floor((H)^(1/3))
            xp4bound = floor((H)^(1/3))
            L.append((H,(x2,x4bound,xp4bound,x6),[0,3]))

        if y4 == yp4 and y4 < min(y2,y6):
            x4 += 1
            xp4 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            x6bound = floor((H)^(1/2))
            L.append((H,(x2bound,x4,xp4,x6bound),[1,2]))

        if y4 == y6 and y4 < min(y2,yp4):
            x4 += 1
            x6 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            xp4bound = floor((H)^(1/3))
            L.append((H,(x2bound,x4,xp4bound,x6),[1,3]))

        if yp4 == y6 and yp4 < min(y2,y4):
            xp4 += 1
            x6 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            x4bound = floor((H)^(1/3))
            L.append((H,(x2bound,x4bound,xp4,x6),[2,3]))

        if y2 == y4 and y4 == yp4 and y2 < y6:
            x2 += 1
            x4 += 1
            xp4 += 1
            H = y2
            x6bound = floor((H)^(1/2))
            L.append((H,(x2,x4,xp4,x6bound),[0,1,2]))

        if y2 == y4 and y4 == y6 and y2 < yp4:
            x2 += 1
            x4 += 1
            x6 += 1
            H = y2
            xp4bound = floor((H)^(1/3))
            L.append((H,(x2,x4,xp4bound,x6),[0,1,3]))

        if y2 == yp4 and y2 == y6 and y2 < y4:
            x2 += 1
            xp4 += 1
            x6 += 1
            H = y2
            x4bound = floor((H)^(1/3))
            L.append((H,(x2,x4bound,x4,x6),[0,2,3]))

        if y4 == yp4 and y4 == y6 and y4 < y2:
            x4 += 1
            xp4 += 1
            x6 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            L.append((H,(x2bound,x4,xp4,x6),[1,2,3]))

        if y2 == y4 and y2 == yp4 and y2 == y6:
            x2 += 1
            x4 += 1
            xp4 += 1
            x6 += 1
            H = y2
            L.append((H,(x2,x4,xp4,x6),[0,1,2,3]))
    return L

def height_iterator_two_torsion(M,N):
    L = []
    x2 = floor((M)^(1/6))
    x4 = floor((M)^(1/3))
    H = max([x2^6,x4^3])
    while H <= N:
        y2 = (x2+1)^6
        y4 = (x4+1)^3
        if min([y2,y4]) > N:
            break
        if y2 < y4:
            x2 += 1
            H = y2
            x4bound = floor((H)^(1/3))
            L.append((H,[x2,x4bound],[0]))
        if y4 < y2:
            x4 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            L.append((H,[x2bound,x4],[1]))
        if y2 == y4:
            x2 += 1
            x4 += 1
            H = y2
            L.append((H,[x2,x4],[0,1]))
    return L

def height_iterator_three_torsion(M,N):
    L = []
    x1 = floor((M)^(1/12))
    x3 = floor((M)^(1/4))
    H = max([x1^12,x3^4])
    while H <= N:
        y1 = (x1+1)^12
        y3 = (x3+1)^4
        if min([y1,y3]) > N:
            break
        if y1 < y3:
            x1 += 1
            H = y1
            x3bound = floor((H)^(1/4))
            L.append((H,[x1,x3bound],[0]))
        if y3 < y1:
            x3 += 1
            H = y3
            x1bound = floor((H)^(1/12))
            L.append((H,[x1bound,x3],[1]))
        if y1 == y3:
            x1 += 1
            x3 += 1
            H = y1
            L.append((H,[x1,x3],[0,1]))
    return L

def height_iterator_full_weierstrass(M,N):
    raise NotImplementedError("Not yet implemented.")


def coefficients_from_height(H,coeffs,indices,model):
    L = []
    for S in list(powerset(indices))[1:]:
        B = []
        for j in range(len(coeffs)):
	    if j in S:
                B.append([-coeffs[j],coeffs[j]])
            elif j in indices:
                B.append(srange(-coeffs[j]+1,coeffs[j]))
	    else:
                B.append(srange(-coeffs[j],coeffs[j]+1))
        C = CartesianProduct(*B).list()
        for c in C:
            L.append(c)

    L2 = []
    for c in L:
        C = coeffs_to_a_invariants(c,model)
        if not is_singular(C):
            L2.append((H,C))
    return L2

def coeffs_to_a_invariants(c,model):
    if model == "short_weierstrass":
        C = [0,0,0,c[0],c[1]]
    elif model == "rank_one":
        C = [0,c[0],c[1],c[2],0]
    elif model == "rank_two":
        a2,a4,b4,a6 = c[0],c[1],c[2],c[3]
        I = 3*a4^2+b4^2-3*a2*a6
        J = -27/4*a2^2*a4^2+18*a4^2*b4 - 2*b4^3 + 9*a2*b4*a6 - 27*a6^2
        if J.denominator()==4:
            I = I*16
            J = J*64
        C = [0,0,0,-27*I,-27*J]
    elif model == "two_torsion":
        C = [0,c[0],0,c[1],0]
    elif model == "three_torsion":
        C = [c[0],0,c[1],0,0]
    elif model == "full_weierstrass":
        raise NotImplementedError("Not yet implemented.")
    else: raise IOError("Please enter recognized Weierstrass family of curves.")

    return C

def is_singular(C):
    a1 = C[0]
    a2 = C[1]
    a3 = C[2]
    a4 = C[3]
    a6 = C[4]

    b2 = a1^2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3^2 + 4*a6
    b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2

    Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
    return Delta==0
    

def set_magma_class_group_bounds(proof = True):
    if proof == False:
        magma.eval('SetClassGroupBounds("GRH")')
    elif proof == True:
        magma.eval('SetClassGroupBounds("PARI")')


def data_by_height(cbh, invariant = 'two_selmer', proof = True):
#INPUT: cbh = List of (height,(a1,a2,a3,a4,a6))
#       invariant = 'rank', 'two_selmer_rank', 'two_selmer_size', 'two_torison_rank', 'reduced_two_selmer_rank', 'three_selmer_rank', 'three_selmer_size', 'three_torsion_rank', or 'reduced_three_selmer_rank'
#OUTPUT: List of (height,(a1,a2,a3,a4,a6),invariant)
    
    t=time.time()
    output = []
    problems = []
    if invariant == 'rank':
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],E.rank(use_database=True, only_use_mwrank = False)))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems

    elif invariant == 'sha_order':
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],E.sha().an_numerical(prec=14,proof=False)))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    elif invariant == 'two_selmer_rank':
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],E.selmer_rank()))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    elif invariant == 'two_selmer_size':
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],2^(E.selmer_rank())))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    elif invariant == 'two_torsion_rank':
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],E.two_torsion_rank()))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems

    
    elif invariant == 'reduced_two_selmer_rank':
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],E.selmer_rank()-E.two_torsion_rank()))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    elif invariant == 'reduced_two_selmer_size':
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],2^(E.selmer_rank()-E.two_torsion_rank())))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    elif invariant == 'three_selmer_rank':
        set_magma_class_group_bounds(proof)
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],E.three_selmer_rank()))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems

    elif invariant == 'three_selmer_size':
        set_magma_class_group_bounds(proof)
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                output.append((C[0],C[1],3^(E.three_selmer_rank())))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    #Below we assume the ground field doesn't contain cube roots of unity.
    
    elif invariant == 'three_torsion_rank':
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                ntors = E.torsion_order()
                if ntors % 3 == 0:
                    output.append((C[0],C[1],1))
                else:
                    output.append((C[0],C[1],0))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    elif invariant == 'reduced_three_selmer_rank':
        set_magma_class_group_bounds(proof)
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                ntors = E.torsion_order()
                if ntors % 3 == 0:
                    output.append((C[0],C[1],E.three_selmer_rank() - 1))
                else:
                    output.append((C[0],C[1],E.three_selmer_rank()))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    elif invariant == 'reduced_three_selmer_rank':
        set_magma_class_group_bounds(proof)
        for C in cbh:
            try:
                E = EllipticCurve(C[1])
                ntors = E.torsion_order()
                if ntors % 3 == 0:
                    output.append((C[0],C[1],3^(E.three_selmer_rank() - 1)))
                else:
                    output.append((C[0],C[1],3^(E.three_selmer_rank())))
            except:
                problems.append(C)
        print(time.time()-t)
        return output,problems
    
    
    else:
        raise NotImplementedError('Invariant not yet Implemented')

def good_primes_by_height(cbh,start,stop):
    t=time.time()
    output = []
    problems = []
    for C in cbh:
        E = EllipticCurve(C[1])
        L = []
	for p in primes(start,stop):
            S = [ZZ(E.is_ordinary(p)),
                 ZZ(E.has_good_reduction(p)),
                 ZZ(E.galois_representation().is_surjective(p))]
            L.append(S)
        output.append((C[0],C[1],L))
    print(time.time()-t)
    return output

def averaged_data(L, filename, return_data=False):
    """
    INPUT:
     -- L: list of tuples (H, coeffs, data), where
        H: height
        coeffs: list of coefficients of elliptic curve giving that height
        data: the type of data that's being averaged, e.g., 2-Selmer, 3-Selmer
     -- filename: name of output file, passed in as a string

    OUTPUT:
     -- writes the averaged data to file
     -- (optional) returns numpy array of dimension len(list) x 2,
        where the first column has heights,
        and the second the averaged data up to that height
    """
    X = np.array([C[0] for C in L])
    V = np.array([C[2] for C in L])
    N = np.arange(1,X.shape[0]+1,dtype=np.float64)
    Y = np.cumsum(V)/N
    I = X[:-1] != X[1:]
    I = np.append(I,True)
    Z = np.vstack([X[I],Y[I]]).T
    np.savetxt(filename, Z)
    if return_data:
       return Z

@parallel
def compute_data(height_bound,model,invariant,filename1="output.txt",\
                 filename2="raw_data.txt",filename3="problems.txt",
		 proof=True):
    L1 = height_iterator(0,height_bound,model)
    L2 = []
    for C in L1:
        L2 += coefficients_from_height(C[0],C[1],C[2],model)

    L3 = data_by_height(L2,invariant,proof)
    
    averaged_data(L3[0],filename1)

    L4 = [flatten(C) for C in L3[0]]
    np.savetxt(filename2,L4)
    L5 = [flatten(C) for C in L3[1]]
    np.savetxt(filename3,L5)