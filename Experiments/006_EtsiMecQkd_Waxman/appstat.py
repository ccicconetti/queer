#!/usr/bin/env python3

import sys

sum_weight = 0.0
sum_rate = 0.0
sum_load = 0.0
for line in sys.stdin:
    tokens = line.rstrip().split(",")
    assert len(tokens) == 4
    (weight, load, rate) = (float(tokens[1]), float(tokens[2]), float(tokens[3]))
    sum_weight += weight
    sum_rate += weight * rate
    sum_load += weight * load

print(f"avg rate: {sum_rate / sum_weight}")
print(f"avg load: {sum_load / sum_weight}")
