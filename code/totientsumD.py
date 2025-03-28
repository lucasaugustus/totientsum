#! /usr/bin/env python3

from utils import *
from time import time
from math import log2, log10

def totientsum_D(n):
    if n <= 10: return 0 if n < 0 else (0,1,2,4,6,10,12,18,22,28,32)[n]
    
    a = introot(n**2, 3)
    b = n // a
    nr = isqrt(n)
    M     = [0] * (   nr + 1)  # M[x]        will store Mertens(x) for small x.
    Mover = [0] * (n//nr + 1)  # Mover[n//x] will store Mertens(x) for large x.
    
    mert = 0
    X, Y, Z = 0, 0, 0
    
    s = nr - (nr == n//nr)
    chi = n // s
    
    
    
    
    for (x, mu) in enumerate(mobiussieve(a+1), start=1):
        v = n // x
        mert += mu
        X += mu * (v * (v+1) // 2)
        
        if x <= nr:
            M[x] = mert
            
            if x > 1:
                for y in range(1, min(b, n//x**2) + 1):
                    Mover[y] -= mu * (n // (y*x))
            
            
            
            
            
        
        elif x == chi:
            if v != b: Mover[v] = mert
            s -= 1
            chi = n // s
        
        if x == a: Z = mert * (b * (b+1) // 2)
    
    for y in range(b, 0, -1):
        v = n // y
        vr = isqrt(v)
        Mv = 1 - v + M[vr] * vr
        for x in range(2, vr+1):
            vx = v // x
            
            Mv -= M[vx] if vx <= nr else Mover[n//vx]
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
    
    methods = (totientsum, \
               #totientsum_C, \
               totientsum_D, \
              )
    
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
            A = m(n)
            print("\t", A)
            answers.append(A)
            print()
            print()
        
        for a in answers: print(a)
        assert len(set(answers)) == 1
    
    if len(argv) > 1:
        print(totientsum_D(int(argv[1])))


