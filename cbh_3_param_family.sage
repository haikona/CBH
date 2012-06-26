import time
import numpy as np

def height_rank_one(M,N):
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
            L.append((H,'x2',[x2,x3bound,x4bound]))

        if y3 < min(y2,y4):
            x3 += 1
            H = y3
            x2bound = floor((H)^(1/6))
            x4bound = floor((H)^(1/3))
            L.append((H,'x3',[x2bound,x3,x4bound]))

        if y4 < min(y2,y3):
            x4 += 1
            H = y4
            x2bound = floor((H)^(1/6))
            x3bound = floor((H)^(1/4))
            L.append((H,'x4',[x2bound,x3bound,x4]))
            
        if y2 == y3 and y2 < y4:
            x2 += 1
            x3 += 1
            H = y2
            x4bound = floor((H)^(1/3))
            L.append((H,'x2x3',[x2,x3,x4bound]))
            
        if y2 == y4 and y2 < y3:
            x2 += 1
            x4 += 1
            H = y2
            x3bound = floor((H)^(1/4))
            L.append((H,'x2x3',[x2,x3bound,x4]))
           
        if y3 == y4 and y3 < y2:
            x3 += 1
            x4 += 1
            H = y2
            x2bound = floor((H)^(1/6))
            L.append((H,'x3x4',[x2bound,x2,x4]))

        if y2 == y3 and y2 == y4:
            x2 += 1
            x3 += 1
            x4 += 1
            H = y2
            L.append((H,'x2x3x4',[x2,x3,x4]))
    return L

def cfh(L):
    L2 = []
    for C in L:
                
        if C[1] == 'x2':
            H = C[2][0]^6
            for j in srange(-C[2][1],C[2][1]+1):
                for k in srange(-C[2][2],C[2][2]+1):
                    L2.append(((-C[2][0],j,k),H))
            for j in srange(-C[2][1],C[2][1]+1):
                for k in srange(-C[2][2],C[2][2]+1):
                    L2.append(((C[2][0],j,k),H))
        
        if C[1] == 'x3':
            H = C[2][1]^4
            for j in srange(-C[2][0],C[2][0]+1):
                for k in srange(-C[2][2],C[2][2]+1):
                    L2.append(((j,-C[2][1],k),H))
            for j in srange(-C[2][0],C[2][0]+1):
                for k in srange(-C[2][2],C[2][2]+1):
                    L2.append(((j,C[2][1],k),H))

        if C[1] == 'x4':
            H = C[2][2]^3
            for j in srange(-C[2][0],C[2][0]+1):
                for k in srange(-C[2][1],C[2][1]+1):
                    L2.append(((j,k,-C[2][2]),H))
            for j in srange(-C[2][0],C[2][0]+1):
                for k in srange(-C[2][1],C[2][1]+1):
                    L2.append(((j,k,C[2][2]),H))
                    
        if C[1] == 'x2x3':
            H = C[2][0]^6
            for j in srange(-C[2][2]+1,C[2][2]):
                L2.append(((-C[2][0],-C[2][1],j),H))
            for j in srange(-C[2][2]+1,C[2][2]):
                L2.append(((-C[2][0],C[2][1],j),H))
            for j in srange(-C[2][2]+1,C[2][2]):
                L2.append(((C[2][0],-C[2][1],j),H))
            for j in srange(-C[2][2]+1,C[2][2]):
                L2.append(((C[2][0],C[2][1],j),H))
                
        if C[1] == 'x2x4':
            H = C[2][0]^6
            for j in srange(-C[2][1]+1,C[2][1]):
                L2.append(((-C[2][0],j,-C[2][2]),H))
            for j in srange(-C[2][1]+1,C[2][1]):
                L2.append(((-C[2][0],j,C[2][2]),H))
            for j in srange(-C[2][2]+1,C[2][2]):
                L2.append(((C[2][0],j,-C[2][2]),H))
            for j in srange(-C[2][2]+1,C[2][2]):
                L2.append(((C[2][0],j,C[2][1]),H))
                
        if C[1] == 'x3x4':
            H = C[2][2]^3
            for j in srange(-C[2][0]+1,C[2][0]):
                L2.append(((j,-C[2][1],-C[2][2]),H))
            for j in srange(-C[2][0]+1,C[2][0]):
                L2.append(((j,-C[2][1],C[2][2]),H))
            for j in srange(-C[2][0]+1,C[2][0]):
                L2.append(((j,C[2][1],-C[2][2]),H))
            for j in srange(-C[2][0]+1,C[2][0]):
                L2.append(((j,C[2][1],C[2][2]),H))
       
        if C[1] =='x2x3x4':
            H = C[2][0]^6
            L2.append(((-C[2][0],-C[2][1],-C[2][2]),H))
            L2.append(((-C[2][0],-C[2][1],C[2][2]),H))
            L2.append(((-C[2][0],C[2][1],-C[2][2]),H))
            L2.append(((C[2][0],-C[2][1],-C[2][2]),H))
            L2.append(((C[2][0],C[2][1],-C[2][2]),H))
            L2.append(((C[2][0],-C[2][1],C[2][2]),H))
            L2.append(((-C[2][0],C[2][1],C[2][2]),H))
            L2.append(((-C[2][0],-C[2][1],-C[2][2]),H))
            
    return L2

@parallel(24)
def avg_h1_selmer(M,N):
    t = time.time()
    L = height_rank_one(M,N)
    L2 = cfh(L)
    
    output = []
    problems = []
    for C in L2:
        try:
            E = EllipticCurve([0,C[0][0],C[0][1],C[0][2],0])
            output.append((E,C[1],E.selmer_rank()))
        except:
            problems.append((C[0]))
    print('Average Selmer size is ', (sum([2^c[2] for c in output])/len(output)).n())
    print(time.time()-t)
    return(output,problems)

@parallel(24)
def average_h1_rank(M,N):
    t = time.time()
    L = height_rank_one(M,N)
    L2 = cfh(L)
    
    output = []
    problems = []
    for C in L2:
        try:
            E = EllipticCurve([0,C[0][0],C[0][1],C[0][2],0])
            output.append((E,C[1],E.rank(use_database=True,only_use_mwrank=False)))
        except:
            problems.append((C[0]))
    print('Average Rank ', (sum([c[2] for c in output])/len(output)).n())
    print(time.time() - t)
    return(output,problems)

@parallel
def h1_rank_list(N):
    L1,L2 = average_h1_rank(0,N)
    X = list_to_array(L1)
    np.savetxt('h1_rank_list_1m.txt',X)
    Y = np.array(L2)
    np.savetxt('h1_rank_problems.txt',Y)

@parallel
def h1_selmer_list(N):
    L1,L2 = avg_h1_selmer(0,N)
    X = list_to_array(L1)
    np.savetxt('h1_selmer_list_1m.txt',X)
    Y = np.array(L2)
    np.savetxt('h1_selmer_problems.txt',Y)