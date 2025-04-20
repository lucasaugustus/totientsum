#! /usr/bin/env python3

from utils import *
from time import time
from math import log2, log10
from itertools import count

def totientsum_J(n):
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
                if len(Mblock) == 0: A = v
                Mblock[v] = mert
                B = v
            s -= 1
            chi = n // s
        
        if x == a: Z = mert * (b * (b+1) // 2)
        
        if (x == a and len(Mblock) > 0) or (x > nr and len(Mblock) == b):
            for y in range(1, b+1):
                tmin = 2
                tmax = isqrt(n//y)
                for t in range(tmin, tmax+1):
                    nty = n // (t*y)
                    if nr < nty:
                        if B <= n // nty <= A:
                            Mover[y] -= Mblock[t*y]
                    else: break
            Mblock = {}
    
    for y in range(b, 0, -1):
        v = n // y
        vr = isqrt(v)
        Mv = 0
        for x in count(2):
            if n / (b+1) >= v // x: break
            Mv -= Mover[n//(v//x)]
        # Mv is now Mertens(v).
        Mover[y] += Mv
        Y += y * Mover[y]
    
    return X + Y - Z




















if __name__ == "__main__":
    from sys import argv
    from random import randrange
    from itertools import chain
    from labmath3 import totientsum, totient
    from totientsumC import totientsum_C
    from totientsumD import totientsum_D
    from totientsumE import totientsum_E
    from totientsumF import totientsum_F
    from totientsumG import totientsum_G
    from totientsumH import totientsum_H
    from totientsumI import totientsum_I
    
    methods = (totientsum, \
               #totientsum_C, \
               #totientsum_D, \
               #totientsum_E, \
               #totientsum_F, \     # slow
               #totientsum_G, \
               #totientsum_H, \     # slow
               totientsum_I, \
               totientsum_J, \
              )
    
    if len(argv) >= 2 and argv[1].isdecimal():
        n = int(argv[1])
        print(n)
        m = methods[-1]
        print(m.__name__)
        print(m(n))
        exit()
    
    if "testlow" in argv:
        DEBUG = False
        verbose = False
        real_s = 0
        for n in range(1, int(argv[2])):
            t = totient(n)
            real_s += t
            test_s = methods[-1](n)
            print('\b'*42, n, real_s, test_s, end='', flush=True)
            assert real_s == test_s
        print()
        exit()
    
    verbose = True
    
    if "DEBUG" in argv:
        args = [x for x in argv if x != "DEBUG"]
        DEBUG = True
    else:
        args = argv
        DEBUG = False
    
    numbers = (2**30, 2**33, 2**34, 2**33 * 3, 10**10, 10**11)
    randos = [randrange(10**8, 10**10) for _ in range(5)]
    
    for n in chain(randos, numbers) if len(args) == 1 else [int(args[1])]:
        print()
        print()
        print()
        print()
        print("N ==", n)
        print()
        print()
        answers = []
        for m in methods:
            print(m.__name__)
            print()
            checkpoint = time()
            A = m(n)
            print("\t", A)
            print(time() - checkpoint)
            answers.append(A)
            print()
            print()
        
        for a in answers: print(a)
        assert len(set(answers)) == 1



