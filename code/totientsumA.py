#! /usr/bin/env python3

from utils import *


def totientsieve(limit=inf):
    """
    Uses a segmented sieve to compute the totients strictly less than limit.
    
    The time- and space-complexities to iterate over the first n terms
    are within logarithmic factors of O(n) and O(sqrt(n)), respectively.
    
    Input: limit -- an integer.  Default = inf.
    
    Output: Sequence of integers
    
    Examples:
    
    >>> list(totientsieve(21))
    [1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18, 8]
    """
    if limit <= 1: return
    yield 1
    pg = primegen()
    primes = [next(pg)]
    nextp = next(pg)
    lo, hi = 2, min(nextp**2, limit)
    # We can sieve up to hi - 1.
    while lo < limit:
        ints = list(range(lo, hi))
        tots = [n for n in range(lo, hi)]
        for p in primes:
            for n in range((-lo) % p, hi - lo, p):
                ints[n] //= p
                tots[n] -= tots[n] // p
            pp = p*p
            while pp < hi:
                for n in range((-lo) % pp, hi - lo, pp):
                    ints[n] //= p
                pp *= p
        # Any entries in ints that are not 1 are prime divisors of their
        # corresponding numbers that were too large to be sieved out.
        for n in range(hi - lo):
            p = ints[n]
            if p != 1:
                tots[n] //= p
                tots[n] *= p - 1
        
        yield from tots
        
        primes.append(nextp)
        nextp = next(pg)
        lo, hi = hi, min(nextp**2, limit)


def totientsum_A(N): return sum(totientsieve(N+1))


if __name__ == "__main__":
    from sys import argv
    print(totientsum_A(int(argv[1])))


