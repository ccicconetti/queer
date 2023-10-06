#!/usr/bin/env python3

import sys

sum_weight = 0.0
sum_rate = 0.0
sum_load = 0.0
n = 0 
for line in sys.stdin:
    tokens = line.rstrip().split(",")
    assert len(tokens) == 4
    (weight, load, rate) = (float(tokens[1]), float(tokens[2]), float(tokens[3]))
    sum_weight += weight
    sum_rate += rate
    sum_load += load
    n += 1

print(f"avg weight: {sum_weight / n}")
print(f"avg rate:   {sum_rate / n}")
print(f"avg load:   {sum_load / n}")
