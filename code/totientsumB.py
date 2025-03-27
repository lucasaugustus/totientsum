#! /usr/bin/env python3

from utils import *


def totientsum_B(n):
    """
    This algorithm is derived by using the Dirichlet hyperbola method on phi = Id * mu.
    The time complexity is O(n^(3/4))-ish: This is due to our use of the Deleglise-Rivat algorithm for computing the
    Mertens function; its O(n^(2/3)) time-complexity makes the optimal choice of parameter a O(n^(3/4)).
    The overall space complexity is O(n^(3/8))-ish, due to the segmented Mobius sieve in stage X.
    Space usage from the Mertens calls is O(n^(1/3))-ish.
    
    Using the Helfgott-Thompson algorithm for the Mertens calls would bring the time complexity down to O(n^(5/7)),
    and the space complexity would be brought down to O(n^(5/14)) == O(n^(0.3571428...)).
    The space complexity from the Mertens calls would then be O(n^0.3).
    
    Further space savings could be obtained by using the Helfgott sieve for the Mobius function:
    the sieving stage would use O(N^(1/4)) space, and the overall space complexity would be controlled by the Mertens
    calls (O(N^(1/3)) if using Deleglise-Rivat, or O(N^(3/10)) if using Helfgott-Thompson).
    """
    if n < 3: return 0 if n < 0 else (0, 1, 2)[n]
    a = introot(n**3, 4)
    b = n // a
    Z = mertens(a) * (b * (b+1) // 2)
    Y = sum(y * mertens(n//y) for y in range(1, b+1))
    X = sum(mu * ((n//x) * ((n//x)+1) // 2) for (x, mu) in enumerate(mobiussieve(a+1), start=1))
    return X + Y - Z


if __name__ == "__main__":
    from sys import argv
    print(totientsum_B(int(argv[1])))


