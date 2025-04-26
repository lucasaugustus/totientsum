#! /usr/bin/env python3

from utils import *


def totientsumM_A(N): return sum(totientsieve(N+1))


if __name__ == "__main__":
    from sys import argv
    print(totientsumM_A(int(argv[1])))


