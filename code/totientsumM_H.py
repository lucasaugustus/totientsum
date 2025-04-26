#! /usr/bin/env python3

from utils import *
from time import time
from math import log2, log10

def totientsumM_H(n):
    if n <= 10: return 0 if n < 0 else (0,1,2,4,6,10,12,18,22,28,32)[n]
    
    a = introot(n**2, 3)
    b = n // a
    nr = isqrt(n)
    Mover = [0] * (b + 1)  # Mover[n//x] will store Mertens(x) for large x.
    Mblock = {}
    
    mert = 0
    X, Y, Z = 0, 0, 0
    
    s = nr - (nr == n//nr)
    chi = n // s
    
    d = b
    xp = isqrt(n//d)
    
    for (x, mu) in enumerate(mobiussieve(a+1), start=1):
        v = n // x
        mert += mu
        X += mu * (v * (v+1) // 2)
        
        if x <= nr:
            Mblock[x] = mert
            
            if x > 1:
                for y in range(1, min(b, n//x**2) + 1):
                    Mover[y] -= mu * (n // (y*x))
            
            while x == xp:
                Mover[d] += 1 - (n//d) + x * mert
                d = (d - 1) if d > 1 else n
                xp = isqrt(n//d)
            
            if x % b == 0 or x == nr:
                A = 1 + (b * (x//b) if x % b != 0 else (x - b))
                assert A == min(Mblock)
                for t in range(1, b+1):
                    lmin = 1 + n // (t * (x+1))
                    lmax = min(isqrt(n//t), n // (t * A))
                    for l in range(lmin, lmax + 1):
                        k = n // (l*t)
                        assert A <= k <= x
                        Mover[t] -= Mblock[k]
                Mblock = {}
        
        elif x == chi:
            if v != b:
                for y in range(1, b+1):
                    tmin = 2
                    tmax = isqrt(n//y)
                    for t in range(tmin, tmax+1):
                        if nr < n//(t*y):
                            if n//x == n // (n//(t*y)):
                                Mover[y] -= mert
            s -= 1
            chi = n // s
        
        if x == a: Z = mert * (b * (b+1) // 2)
    
    for y in range(b, 0, -1):
        v = n // y
        vr = isqrt(v)
        Mv = 0
        for x in range(2, vr+1):
            vx = v // x
            if vx > nr:
                if n//vx <= b:
                    Mv -= Mover[n//vx]
        # Mv is now Mertens(v).
        Mover[y] += Mv
        Y += y * Mover[y]
    
    return X + Y - Z




if __name__ == "__main__":
    from sys import argv
    print(totientsumM_H(int(argv[1])))


