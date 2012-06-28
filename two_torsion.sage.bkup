import time
import numpy as np

def height_two_torsion(M,N):
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
            L.append((H,'x2',x2,x4bound))
        if y4 < y2:
            x4 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            L.append((H,'x4',x2bound,x4))
        if y2 == y4:
            x2 += 1
            x4 += 1
            H = y2
            L.append((H,'x2x4',x2,x4))
    return L
    
    
def coeffs_from_height2(L):
    L2 = []
    for C in L:
        
        if C[1] == 'x2':
            H = C[2]^6
            for j in srange(-C[3],C[3]+1):
                L2.append(((-C[2],j),H))
            for j in srange(-C[3],C[3]+1):
                L2.append(((C[2],j),H))
        
        if C[1] == 'x4':
            H = C[3]^3
            for k in srange(-C[2],C[2]+1):
                L2.append(((k,-C[3]),H))
            for k in srange(-C[2],C[2]+1):
                L2.append(((k,C[3]),H))

        if C[1] == 'x2x4':
            H = C[2]^6
            for j in srange(-C[3]+1,C[3]):
                L2.append(((-C[2],j),H))
            for j in srange(-C[3]+1,C[3]):
                L2.append(((C[2],j),H))
            for k in srange(-C[2]+1,C[2]):
                L2.append(((k,-C[3]),H))
            for k in srange(-C[2]+1,C[2]):
                L2.append(((k,C[3]),H))
            L2.append(((-C[2],-C[3]),H))
            L2.append(((-C[2],C[3]),H)) 
            L2.append(((C[2],-C[3]),H))
            L2.append(((C[2],C[3]),H))   
    return L2


def avg_htwotor_sel2(M,N):
    t = time.time()
    L = height_two_torsion(M,N)
    L2 = coeffs_from_height2(L)
    
    output = []
    problems = []
    for C in L2:
        try:
            E = EllipticCurve([0,C[0][0],0,C[0][1],0])
            output.append((E,C[1],E.selmer_rank()))
        except:
            problems.append((C[0]))
    print(time.time()-t)
    print('Average 2-Selmer size is ', (sum([2^c[2] for c in output])/len(output)).n())
    return(output,problems)



def average_htwotor_rank(M,N):
    t = time.time()
    L = height_two_torsion(M,N)
    L2 = coeffs_from_height2(L)
    
    output = []
    problems = []
    for C in L2:
        try:
            E = EllipticCurve([0,C[0][0],0,C[0][1],0])
            output.append((E,C[1],E.rank(use_database=True,only_use_mwrank=False)))
        except:
            problems.append((C[0]))
    print(time.time()-t)
    print('Average Rank ', (sum([c[2] for c in output])/len(output)).n())
    return(output,problems)



def htwotor_rank_list(N):
    L1,L2 = average_htwotor_rank(0,N)
    X = list_to_array(L1)
    np.savetxt('htwotor_rank_list_1m.txt',X)
    Y = np.array(L2)
    np.savetxt('htwotor_rank_problems.txt',Y)


def htwotor_sel2_list(N):
    L1,L2 = avg_htwotor_sel2(0,N)
    X = list_to_array(L1)
    np.savetxt('htwotor_sel2_list_1m.txt',X)
    Y = np.array(L2)
    np.savetxt('htwotor_sel2_problems.txt',Y)