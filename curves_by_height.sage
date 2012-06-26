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
    raise NotImplementedError("Not yet implemented.")

def height_iterator_two_torsion(M,N):
    raise NotImplementedError("Not yet implemented.")

def height_iterator_three_torsion(M,N):
    raise NotImplementedError("Not yet implemented.")

def height_iterator_full_weierstrass(M,N):
    raise NotImplementedError("Not yet implemented.")