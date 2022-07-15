#!/usr/bin/env python3

import math


def l_i(F_i: float, F: float):
    return math.floor(math.log((4 * F_i - 1) / 3) / math.log((4 * F - 1) / 3))


def purif(x1: float, x2: float):
    return x1 * x2 / (x1 * x2 + (1 - x1) * (1 - x2))


def purif_round(x: float, n: int):
    assert n >= 1
    ret = x
    for r in range(n):
        ret = purif(x, ret)
    return ret


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

    with open("fidelity_round.dat", "w") as outfile:
        for f in range(50, 101, 1):
            values = []
            for r in range(1, 5, 1):
                values.append(f"{purif_round(f / 100.0, r)}")
            outfile.write(f'{f/100.0} {" ".join(values)}\n')
