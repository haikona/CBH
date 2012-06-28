import time
import numpy as np

def height_iterator(M,N):
    L = []

    a = floor((M/4)^(1/3))
    b = floor((M/27)^(1/2))
    H = max([4*a^3,27*b^2])
    while H <= N:
        a2 = 4*(a+1)^3
        b2 = 27*(b+1)^2
        if min([a2,b2]) > N:
            break
        if a2 < b2:
            a += 1
            H = a2
            bbound = floor((H/27)^(1/2))
            L.append((H,'a',a,bbound))
        elif b2 < a2:
            b += 1
            H = b2
            abound = floor((H/4)^(1/3))
            L.append((H,'b',a,b))
        elif a2 == b2:
            a += 1
            b += 1
            H = a2
            L.append((H,'ab',a,b))
    return L
    
def coeffs_from_height(L):
    L2 = []
    for C in L:
        
        if C[1] == 'a':
            H = 4*C[2]^3
            for j in srange(-C[3],C[3]+1):
                L2.append(((-C[2],j),H))
            for j in srange(-C[3],C[3]+1):
                L2.append(((C[2],j),H))
        
        if C[1] == 'b':
            H = 27*C[3]^2
            for k in srange(-C[2],C[2]+1):
                L2.append(((k,-C[3]),H))
            for k in srange(-C[2],C[2]+1):
                L2.append(((k,C[3]),H))

        if C[1] == 'ab':
            H = 4*C[2]^3
            for j in srange(-C[3]+1,C[3]):
                L2.append(((-C[2],j),H))
            for j in srange(-C[3]+1,C[3]):
                L2.append(((C[2],j),H))
            for k in srange(-C[2]+1,C[2]):
                L2.append(((k,-C[3]),H))
            for k in srange(-C[2]+1,C[2]):
                L2.append(((k,C[3]),H))
            L2.append(((C[2],-C[3]),H))
            L2.append(((C[2],C[3]),H))  
    return L2

def two_selmer_size_by_height(M,N):
    t = time.time()
    L = height_iterator(M,N)
    L2 = coeffs_from_height(L)

    output = []
    problems = []
    for C in L2:
        E = EllipticCurve(C[0])
        try:
            output.append((E,C[1],2^(E.selmer_rank())))
        except:
            problems.append((E,C[1]))
    print(time.time()-t)
    return (output,problems)

def ranks_by_height(M,N):
    t = time.time()
    L = height_iterator(M,N)
    L2 = coeffs_from_height(L)
    
    output = []
    problems = []
    for C in L2:
        E = EllipticCurve(C[0])
        try:
            output.append((E,C[1],E.rank(use_database=True,only_use_mwrank=False)))
        except:
            problems.append((E,C[1]))
    print(time.time()-t)
    return (output,problems)

def list_to_array(L):
    X = np.array([[C[-2],C[-1]] for C in L])
    return X

def list_to_array_with_coeffs(L):
    X = np.array([[C[0].a4(),C[0].a6(),C[1],C[2]] for C in L])
    return X

def average_ranks(L):
    X = L[:,-2]
    N = np.arange(1,X.shape[0]+1,dtype=np.float64)
    Y = np.cumsum(L[:,-1])/N

    I = X[:-1] != X[1:]
    I = np.append(I,True)

    Z = np.vstack([X[I],Y[I]]).T
    return Z

def rank_list(N):
    L1,L2 = ranks_by_height(0,N)
    X = list_to_array_with_coeffs(L1)
    np.savetxt('rank_list_1m.txt',X)
    Y = np.array([[C[0].a4(),C[0].a6(),C[1]] for C in L2])
    np.savetxt('rank_problems.txt',Y)

def selmer_list(N):
    L1,L2 = two_selmer_size_by_height(0,N)
    X = list_to_array_with_coeffs(L1)
    np.savetxt('selmer_list_1m.txt',X)
    Y = np.array([[C[0].a4(),C[0].a6(),C[1]] for C in L2])
    np.savetxt('selmer_problems.txt',Y)

def crunch_ranks(N):
    L1,L2 = ranks_by_height(0,N)
    X = list_to_array(L1)
    Z = average_ranks(X)
    np.savetxt('output.txt',Z)
