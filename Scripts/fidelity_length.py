#!/usr/bin/env python3

import math


def l_i(F_i: float, F: float):
    return math.floor(math.log((4 * F_i - 1) / 3) / math.log((4 * F - 1) / 3))

def purif(x1: float, x2: float):
    return x1 * x2 / (x1 * x2 + (1 - x1) * (1 - x2))


if __name__ == "__main__":
    F = 0.99
    for fidelity in range(93, 100, 1):
        f = fidelity / 100.0
        print(f"{f}: {l_i(f, F)}")

    print(purif(0.88, 0.88))
    print(purif(0.88, 0.92))
    print(purif(0.96, 0.98))
    print(purif(0.92, 0.98))
    print(purif(0.93, 0.97))
