#!/usr/bin/env python3

import sys

even = []
odd = []
cnt = 0
for line in sys.stdin:
    line = line.rstrip()
    if line == "":
        continue

    if cnt % 2 == 0:
        even.append(line)
    else:
        odd.append(line)

    cnt += 1

assert len(even) == len(odd)

for (a,b) in zip(even,odd):
    if a != "0,0" and b != "0,0":
        print(f'{a}\n{b}\n')
