#! /usr/bin/env python3

from utils import *
from time import time
from itertools import count




def totientsum(n):
    if n <= 10: return 0 if n < 0 else (0,1,2,4,6,10,12,18,22,28,32)[n]
    
    a = introot(int((n / log(log(n)))**2), 3)
    b = n // a
    nr = isqrt(n)
    if verbose:
        print("        a:", a)
        print("  sqrt(n):", nr)
        print("        b:", b)
        print("  sqrt(a):", isqrt(a))
        print("sqrt(n)/b:", nr//b)
        print("ln(ln(n)):", log(log(n)))
        starttime = time()
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
        #v = int(str(n // x))    # The int(str( ... )) pushes us back down into the 64-bit data types, when applicable.
        v = n // x
        mert += mu
        X += mu * (v * (v+1) // 2)
        
        if x <= nr:
            Mblock.append(mert)
            
            if x > 1 and mu != 0:
                if verbose and (x<10*b or x%b<2): print("\b"*80, " Phase 1:", x, "%f%%" % (100*x/nr), end='    ', flush=True)
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
                if verbose: print("\b"*80, " P1 Batch", x//b, "%f%%" % (100*x/nr), end='    ', flush=True)
                A = 1 + (b * (x//b) if x % b != 0 else (x - b))
                for t in range(1, b+1):
                    #if verbose and x < 20*b: print("\b"*80, "         ", x//b, t, "%f%%" % (100*x/a), end='    ', flush=True)
                    nt = n // t
                    lmin = 1 + n // (t * (x+1))
                    lmax = min(isqrt(n//t), n // (t * A))
                    for l in range(lmin, lmax + 1):
                        k = nt // l
                        assert A <= k <= x
                        Mover[t] -= Mblock[k - A]
                Mblock.clear()
                if verbose and x == nr:
                    phase1time = time()
                    print("\b"*80 + ("Phase 1 took %f seconds.    " % (phase1time - starttime)))
        
        elif x == chi:
            if verbose and len(Mblock) % 100 == 0: print("\b"*80, " Phase 2:", x, "%f%%" % (100*x/a), end='    ', flush=True)
            if v != b:
                if len(Mblock) == 0: A = v
                Mblock.append(mert)
                B = v
            s -= 1
            chi = n // s
        
        if (x == a and len(Mblock) > 0) or (x > nr and len(Mblock) == b):
            if verbose: print("\b"*80, " P2 Batch", x, "%f%%" % (100*x/a), end='    ', flush=True)
            BnBn = B * n // (B + n)
            for y in range(1, b+1):
                ny = n // y
                tmin = max(2, BnBn // y)
                tmax = min(isqrt(ny), (A + 1) // y)
                for t in range(tmin, tmax+1):
                    nty = ny // t
                    if nr < nty:
                        if B <= n // nty <= A:
                            Mover[y] -= Mblock[A - t*y]
                    else: break
            Mblock.clear()
        
    Z = mert * (b * (b+1) // 2)
    
    if verbose:
        phase2time = time()
        print("\b"*80 + ("Phase 2 took %f seconds.    " % (phase2time - phase1time)))
    
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
    
    if verbose:
        phase3time = time()
        print("\b"*80 + ("Phase 3 took %f seconds." % (phase3time - phase2time)))
    
    return X + Y - Z




if __name__ == "__main__":
    from sys import argv
    print(totientsum(int(argv[1])))




