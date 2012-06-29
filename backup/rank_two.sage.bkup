import time
import numpy as np

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




def coefficients_from_height(L):
    L2 = []
    for C in L:
                
        if C[2] == [0]:
            H = C[1][0]^6
            for j in srange(-C[1][1],C[1][1]+1):
                for k in srange(-C[1][2],C[1][2]+1):
                    for l in srange(-C[1][3],C[1][3]+1):
                        L2.append(((-C[1][0],j,k,l),H))
            for j in srange(-C[1][1],C[1][1]+1):
                for k in srange(-C[1][2],C[1][2]+1):
                    for l in srange(-C[1][3],C[1][3]+1):
                        L2.append(((C[1][0],j,k,l),H))

        if C[2] == [1]:
            H = C[1][1]^3
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][2],C[1][2]+1):
                    for l in srange(-C[1][3],C[1][3]+1):
                        L2.append(((j,-C[1][1],k,l),H))
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][2],C[1][2]+1):
                    for l in srange(-C[1][3],C[1][3]+1):
                        L2.append(((j,C[1][1],k,l),H))

        if C[2] == [2]:
            H = C[1][2]^3
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][1],C[1][1]+1):
                    for l in srange(-C[1][3],C[1][3]+1):
                        L2.append(((j,k,-C[1][2],l),H))
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][1],C[1][1]+1):
                    for l in srange(-C[1][3],C[1][3]+1):
                        L2.append(((j,k,C[1][2],l),H))

        if C[2] == [3]:
            H = C[1][3]^2
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][1],C[1][1]+1):
                    for l in srange(-C[1][2],C[1][2]+1):
                        L2.append(((j,k,l,-C[1][3]),H))
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][1],C[1][1]+1):
                    for l in srange(-C[1][2],C[1][2]+1):
                        L2.append(((j,k,l,C[1][3]),H))
                        
                    
        if C[2] == [0,1]:
            H = C[1][0]^6
            for j in srange(-C[1][2],C[1][2]+1):
                for k in srange(-C[1][3],C[1][3]+1):
                    L2.append(((-C[1][0],-C[1][1],j,k),H))
                    L2.append(((-C[1][0],C[1][1],j,k),H))
                    L2.append(((C[1][0],-C[1][1],j,k),H))
                    L2.append(((C[1][0],C[1][1],j,k),H))
                
        if C[2] == [0,2]:
            H = C[1][0]^6
            for j in srange(-C[1][1],C[1][1]+1):
                for k in srange(-C[1][3],C[1][3]+1):
                    L2.append(((-C[1][0],j,-C[1][2],k),H))
                    L2.append(((-C[1][0],j,C[1][2],k),H))
                    L2.append(((C[1][0],j,-C[1][2],k),H))
                    L2.append(((C[1][0],j,C[1][2],k),H))

        if C[2] == [0,3]:
            H = C[1][0]^6
            for j in srange(-C[1][1],C[1][1]+1):
                for k in srange(-C[1][2],C[1][2]+1):
                    L2.append(((-C[1][0],j,k,-C[1][3]),H))
                    L2.append(((-C[1][0],j,k,C[1][3]),H))
                    L2.append(((C[1][0],j,k,-C[1][3]),H))
                    L2.append(((C[1][0],j,k,C[1][3]),H))
                        
        if C[2] == [1,2]:
            H = C[1][1]^3
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][3],C[1][3]+1):
                    L2.append(((j,-C[1][1],-C[1][2],k),H))
                    L2.append(((j,-C[1][1],C[1][2],k),H))
                    L2.append(((j,C[1][1],-C[1][2],k),H))
                    L2.append(((j,C[1][1],C[1][2],k),H))

        if C[2] == [1,3]:
            H = C[1][1]^3
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][2],C[1][2]+1):
                    L2.append(((j,-C[1][1],k,-C[1][3]),H))
                    L2.append(((j,-C[1][1],k,C[1][3]),H))
                    L2.append(((j,C[1][1],k,-C[1][3]),H))
                    L2.append(((j,C[1][1],k,C[1][3]),H))

        if C[2] == [2,3]:
            H = C[1][2]^3
            for j in srange(-C[1][0],C[1][0]+1):
                for k in srange(-C[1][1],C[1][1]+1):
                    L2.append(((j,k,-C[1][2],-C[1][3]),H))
                    L2.append(((j,k,-C[1][2],C[1][3]),H))
                    L2.append(((j,k,C[1][2],-C[1][3]),H))
                    L2.append(((j,k,C[1][2],C[1][3]),H))
                        
        if C[2] == [0,1,2]:
            H = C[1][0]^6
            for j in srange(-C[1][3],C[1][3]+1):
                L2.append(((-C[1][0],-C[1][1],-C[1][2],j),H))
                L2.append(((-C[1][0],-C[1][1],C[1][2],j),H))
                L2.append(((-C[1][0],C[1][1],-C[1][2],j),H))
                L2.append(((C[1][0],-C[1][1],-C[1][2],j),H))
                L2.append(((C[1][0],C[1][1],-C[1][2],j),H))
                L2.append(((C[1][0],-C[1][1],C[1][2],j),H))
                L2.append(((-C[1][0],C[1][1],C[1][2],j),H))
                L2.append(((C[1][0],C[1][1],C[1][2],j),H))
                
        if C[2] == [0,1,3]:
            H = C[1][0]^6
            for j in srange(-C[1][2],C[1][2]+1):
                L2.append(((-C[1][0],-C[1][1],j,-C[1][3]),H))
                L2.append(((-C[1][0],-C[1][1],j,C[1][3]),H))
                L2.append(((-C[1][0],C[1][1],j,-C[1][3]),H))
                L2.append(((C[1][0],-C[1][1],j,-C[1][3]),H))
                L2.append(((C[1][0],C[1][1],j,-C[1][3]),H))
                L2.append(((C[1][0],-C[1][1],j,C[1][3]),H))
                L2.append(((-C[1][0],C[1][1],j,C[1][3]),H))
                L2.append(((C[1][0],C[1][1],j,C[1][3]),H))
                
        if C[2] == [0,2,3]:
            H = C[1][0]^6
            for j in srange(-C[1][1],C[1][1]+1):
                L2.append(((-C[1][0],j,-C[1][2],-C[1][3]),H))
                L2.append(((-C[1][0],j,-C[1][2],C[1][3]),H))
                L2.append(((-C[1][0],j,C[1][2],-C[1][3]),H))
                L2.append(((C[1][0],j,-C[1][2],-C[1][3]),H))
                L2.append(((C[1][0],j,C[1][2],-C[1][3]),H))
                L2.append(((C[1][0],j,-C[1][2],C[1][3]),H))
                L2.append(((-C[1][0],j,C[1][2],C[1][3]),H))
                L2.append(((C[1][0],j,C[1][2],C[1][3]),H))
                                                    
        if C[2] == [1,2,3]:
            H = C[1][1]^3
            for j in srange(-C[1][0],C[1][0]+1):
                L2.append(((j,-C[1][1],-C[1][2],-C[1][3]),H))
                L2.append(((j,-C[1][1],-C[1][2],C[1][3]),H))
                L2.append(((j,-C[1][1],C[1][2],-C[1][3]),H))
                L2.append(((j,C[1][1],-C[1][2],-C[1][3]),H))
                L2.append(((j,C[1][1],C[1][2],-C[1][3]),H))
                L2.append(((j,C[1][1],-C[1][2],C[1][3]),H))
                L2.append(((j,-C[1][1],C[1][2],C[1][3]),H))
                L2.append(((j,C[1][1],C[1][2],C[1][3]),H))
                
        if C[2] == [0,1,2,3]:
            H = C[1][0]^6
            L2.append(((C[1][0],-C[1][1],-C[1][2],-C[1][3]),H))
            L2.append(((C[1][0],-C[1][1],-C[1][2],C[1][3]),H))
            L2.append(((C[1][0],-C[1][1],C[1][2],-C[1][3]),H))
            L2.append(((C[1][0],C[1][1],-C[1][2],-C[1][3]),H))
            L2.append(((C[1][0],C[1][1],C[1][2],-C[1][3]),H))
            L2.append(((C[1][0],C[1][1],-C[1][2],C[1][3]),H))
            L2.append(((C[1][0],-C[1][1],C[1][2],C[1][3]),H))
            L2.append(((C[1][0],C[1][1],C[1][2],C[1][3]),H))    
            L2.append(((-C[1][0],-C[1][1],-C[1][2],-C[1][3]),H))
            L2.append(((-C[1][0],-C[1][1],-C[1][2],C[1][3]),H))
            L2.append(((-C[1][0],-C[1][1],C[1][2],-C[1][3]),H))
            L2.append(((-C[1][0],C[1][1],-C[1][2],-C[1][3]),H))
            L2.append(((-C[1][0],C[1][1],C[1][2],-C[1][3]),H))
            L2.append(((-C[1][0],C[1][1],-C[1][2],C[1][3]),H))
            L2.append(((-C[1][0],-C[1][1],C[1][2],C[1][3]),H))
            L2.append(((-C[1][0],C[1][1],C[1][2],C[1][3]),H)) 

    return L2


def I(a2,a4,b4,a6):
    return(3*a4^2+b4^2-3*a2*a6)


def J(a2,a4,b4,a6):
    return(-27/4*a2^2*a4^2+18*a4^2*b4 - 2*b4^3 + 9*a2*b4*a6 - 27*a6^2)


def EC(v):
    IE = I(v[0],v[1],v[2],v[3])
    JE = J(v[0],v[1],v[2],v[3])
    return(EllipticCurve([-27*IE,-27*JE]))


def avg_sel2(M,N):
    t = time.time()
    L = height_iterator_rank_two(M,N)
    L2 = coefficients_from_height(L)
        
    output = []
    problems = []
    for C in L2:
        try:
            E = EC(C[0])
            E2 = E.integral_model()
            output.append((E2,C[1],E2.selmer_rank()))
        except:
            problems.append((C[0]))
    print(time.time()-t)
    print('Average Selmer size is ', (sum([2^c[2] for c in output])/len(output)).n())
    return(output,problems)