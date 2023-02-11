#!/usr/bin/env python3

import sys

data = []
for line in sys.stdin:
    data.append(float(line.rstrip()))

sum = 0
ssum = 0
for x in data:
    sum += x
    ssum += x * x

print(f"{sum * sum / (len(data) * ssum)}")
