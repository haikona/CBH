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
            L.append((H,'x1',x1,x3bound))
        if y3 < y1:
            x3 += 1
            H = y3
            x1bound = floor((H)^(1/12))
            L.append((H,'x3',x1bound,x3))
        if y1 == y3:
            x1 += 1
            x3 += 1
            H = y1
            L.append((H,'x1x3',x1,x3))
    return L