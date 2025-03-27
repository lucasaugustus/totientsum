#! /usr/bin/env python3

from utils import *
from time import time
from math import log2, log10

def totientsum_C(n):
    if n <= 10: return 0 if n < 0 else (0,1,2,4,6,10,12,18,22,28,32)[n]
    
    a = introot(n**2, 3)
    b = n // a
    nr = isqrt(n)
    M     = [0] * (   nr + 1)  # M[k]        will store Mertens(k) for small k.
    Mover = [0] * (n//nr + 1)  # Mover[x//k] will store Mertens(k) for large k.
    mobs = [0] * (nr+1)
    mert = 0
    X, Y, Z = 0, 0, 0
    
    s = nr - (nr == n//nr)
    Mkey1 = key = n // s
    for (k, mu) in enumerate(mobiussieve(a+1), start=1):
        v = n // k
        mert += mu
        X += mu * (v * (v+1) // 2)
        
        if k <= nr:
            M[k] = mert
            mobs[k] = mu
        
        
        
        
        elif k == key:
            if v != b: Mover[v] = mert
            s -= 1
            key = n // s
        
        if k == a: Z = mert * (b * (b+1) // 2)
    
    #Mover[b] = 0
    
    # Now that we have Mobius values up to xr stored in mobs, and some Mertens values up to y
    # stored in M and Mover, we compute the rest of the needed Mertens values up to x with the formula
    # M(v) == 1 - v - sum(mu(n) * (v//n) + M(v//n)) + isqrt(v) * M(sqrt(v)),
    # where the sum runs over 2 <= n <= sqrt(v).
    
    for y in range(b, 0, -1):
        v = n // y
        vr = isqrt(v)
        Mv = 1 - v + M[vr] * vr
        for k in range(2, vr+1):
            vk = v // k
            Mv -= mobs[k] * vk
            Mv -= M[vk] if vk <= nr else Mover[n//vk]
        # Mv is now Mertens(v).
        Mover[y] += Mv
        Y += y * Mover[y]
    
    # Now we can compute the totient sum.  We use the formula
    # totientsum(n) == 1/2 * sum(mu(k) * (n//k) * (1 + n//k)), where 1 <= k <= n.
    # We exploit the fact that n//k takes many repeated values when sqrt(n) <= k <= n.
    
    return X + Y - Z

















if __name__ == "__main__":
    from sys import argv
    from random import randrange
    from itertools import chain
    from labmath3 import totientsum, totient
    
    
    methods = (totientsum, \
               totientsum_C, \
               
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
        print(totientsum_C(int(argv[1])))


