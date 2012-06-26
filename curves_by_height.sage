# A more systematic attempt at code for generating statistics of invariants
# of elliptic curves ordered by height.

def height_iterator(M,N,model):
    if model == "short_weirstrass":
       	return height_iterator_short_weierstrass(M,N)
    elif model == "rank_one":
        return height_iterator_rank_one(M,N)
    elif model == "two_torsion":
        return height_iterator_two_torsion(M,N)
    elif model == "three_torsion":
        return height_iterator_three_torsion(M,N)
    elif model == "full_weierstrass":
        return height_iterator_full_weierstrass(M,N)
    else: raise InputError("Please enter recognized Weierstrass family of curves.")

def height_iterator_short_weierstrass(M,N):
    raise NotImplementedError("Not yet implemented.")

def height_iterator_rank_one(M,N):
    raise NotImplementedError("Not yet implemented.")

def height_iterator_two_torsion(M,N):
    raise NotImplementedError("Not yet implemented.")

def height_iterator_three_torsion(M,N):
    raise NotImplementedError("Not yet implemented.")

def height_iterator_full_weierstrass(M,N):
    raise NotImplementedError("Not yet implemented.")