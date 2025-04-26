#! /usr/bin/env python3

from utils import *


def totientsumT_A(n):
    if n <= 10: return 0 if n < 0 else (0,1,2,4,6,10,12,18,22,28,32)[n]
    a = introot(n**2, 3)
    b = n // a
    nr = isqrt(n)
    phi = [0] * (nr + 1)
    Phi = [0] * (nr + 1)
    Phiover = [0] * (n//nr + 1)
    X, Y, Z = 0, 0, 0
    
    s = nr - (nr == n // nr)
    key = n // s
    totsum = 0
    for (x, tot) in enumerate(totientsieve(a+1), start=1):
        totsum += tot
        X += tot * (n//x)
        if x <= nr:
            phi[x] = tot
            Phi[x] = totsum
        if x == key:
            if n // x != b: Phiover[n//x] = totsum
            s -= 1
            key = n // s
        if x == a:
            Z = n * (n+1) // 2 + b * totsum
    
    for y in range(b, 1, -1):
        v = n // y
        vr = isqrt(v)
        PhiV = v * (v+1) // 2 + vr * Phi[vr] - v
        for t in range(2, vr + 1):
            q = v // t
            PhiV -= phi[t] * q
            if q <= nr: PhiV -= Phi[q]
            else:       PhiV -= Phiover[n//q]
        Phiover[y] += PhiV
        Y += PhiV
    
    return Z - X - Y




def totientsumT_B(n):
    if n <= 10: return 0 if n < 0 else (0,1,2,4,6,10,12,18,22,28,32)[n]
    a = introot(n**2, 3)
    b = n // a
    nr = isqrt(n)
    phi = [0] * (nr + 1)
    Phi = [0] * (nr + 1)
    Phiover = [0] * (n//nr + 1)
    X, Y, Z = 0, 0, 0
    
    s = nr - (nr == n // nr)
    key = n // s
    totsum = 0
    for (x, tot) in enumerate(totientsieve(a+1), start=1):
        totsum += tot
        X += tot * (n//x)
        if x <= nr:
            phi[x] = tot
            Phi[x] = totsum
        if x == key:
            if n // x != b: Phiover[n//x] = totsum
            s -= 1
            key = n // s
        if x == a:
            Z = n * (n+1) // 2 + b * totsum
    
    for y in range(b, 1, -1):
        v = n // y
        vr = isqrt(v)
        PhiV = v * (v+1) // 2 + vr * Phi[vr] - v
        for t in range(2, vr + 1):
            q = v // t
            PhiV -= phi[t] * q
            if q <= nr: PhiV -= Phi[q]
            else:       PhiV -= Phiover[n//q]
        Phiover[y] += PhiV
        Y += PhiV
    
    return Z - X - Y







if __name__ == "__main__":
    from sys import argv
    from time import time
    from random import randrange
    from itertools import chain
    from labmath3 import totientsum, totient
    
    methods = (totientsum, \
               totientsumT_A, \
               totientsumT_B, \
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



