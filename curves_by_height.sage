# A more systematic attempt at code for generating statistics of invariants
# of elliptic curves ordered by height.

def height_iterator(M,N,model):
    if model == "short_weierstrass":
       	return height_iterator_short_weierstrass(M,N)
    elif model == "rank_one":
        return height_iterator_rank_one(M,N)
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

def height_iterator_two_torsion(M,N):
    raise NotImplementedError("Not yet implemented.")

def height_iterator_three_torsion(M,N):
    raise NotImplementedError("Not yet implemented.")

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
            L.append([H,c])
    
    return L