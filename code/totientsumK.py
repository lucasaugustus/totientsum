#! /usr/bin/env python3

from utils import *
from time import time
from math import log2, log10
from itertools import count

def totientsumM_K(n):
    if n <= 10: return 0 if n < 0 else (0,1,2,4,6,10,12,18,22,28,32)[n]
    
    a = introot(int((n / log(log(n)))**2), 3)
    b = n // a
    nr = isqrt(n)
    if verbose:
        print("      b:", b)
        print("sqrt(n):", nr)
        print("      a:", a)
    Mover = [0] * (b + 1)  # Mover[n//x] will store Mertens(x) for large x.
    Mblock = []
    
    mert = 0
    X, Y, Z = 0, 0, 0
    
    s = nr - (nr == n//nr)
    chi = n // s
    
    d = b
    xp = isqrt(n//d)
    
    for (x, mu) in enumerate(mobiussieve(a+1), start=1):
        #if x == 100: exit()
        v = int(str(n // x))    # The int(str( ... )) pushes us back down into the 64-bit data types, when applicable.
        mert += mu
        X += mu * (v * (v+1) // 2)
        
        if x <= nr:
            Mblock.append(mert)
            
            if x > 1 and mu != 0:
                if verbose:
                    if x % b == 0 or (v//x > b // 100):
                        print("\b"*80, "       ", x, "%f%%" % (100*x/a), end='    ', flush=True)
                if mu > 0:
                    for y in range(1, min(b, v//x) + 1):
                        Mover[y] -= v // y
                else:
                    for y in range(1, min(b, v//x) + 1):
                        Mover[y] += v // y
            
            while x == xp:
                Mover[d] += 1 - (n//d) + x * mert
                d = (d - 1) if d > 1 else n
                xp = isqrt(n//d)
            
            if x % b == 0 or x == nr:
                A = 1 + (b * (x//b) if x % b != 0 else (x - b))
                for t in range(1, b+1):
                    if x // b < 10: print("\b"*80, "       ", x//b, t, "%f%%" % (100*x/a), end='    ', flush=True)
                    nt = n // t
                    lmin = 1 + n // (t * (x+1))
                    lmax = min(isqrt(n//t), n // (t * A))
                    for l in range(lmin, lmax + 1):
                        k = nt // l
                        assert A <= k <= x
                        Mover[t] -= Mblock[k - A]
                Mblock.clear()
        
        elif x == chi:
            if v != b:
                if len(Mblock) == 0: A = v
                Mblock.append(mert)
                B = v
            s -= 1
            chi = n // s
        
        if x == a: Z = mert * (b * (b+1) // 2)
        
        if (x == a and len(Mblock) > 0) or (x > nr and len(Mblock) == b):
            if verbose:
                print("\b"*80, "       ", x, "%f%%" % (100*x/a), end='    ', flush=True)
            for y in range(1, b+1):
                tmin = 2
                tmax = isqrt(n//y)
                for t in range(tmin, tmax+1):
                    nty = n // (t*y)
                    if nr < nty:
                        if B <= n // nty <= A:
                            Mover[y] -= Mblock[A - t*y]
                    else: break
            Mblock.clear()
    
    if verbose: print()
    
    for y in range(b, 0, -1):
        v = n // y
        vr = isqrt(v)
        Mv = 0
        for x in count(2):
            if n >= (b+1) * (v//x): break
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
    from totientsumM_C import totientsumM_C
    from totientsumM_D import totientsumM_D
    from totientsumM_E import totientsumM_E
    from totientsumM_F import totientsumM_F
    from totientsumM_G import totientsumM_G
    from totientsumM_H import totientsumM_H
    from totientsumM_I import totientsumM_I
    from totientsumM_J import totientsumM_J
    
    methods = (totientsum, \
               totientsumM_K, \
              )
    
    verbose = True
    
    if len(argv) >= 2 and argv[1].isdecimal():
        n = int(argv[1])
        print(n)
        m = methods[-1]
        print(m.__name__)
        print(m(n))
        exit()
    
    if "testlow" in argv:
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



