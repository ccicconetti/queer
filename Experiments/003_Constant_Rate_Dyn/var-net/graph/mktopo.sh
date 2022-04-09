#!/bin/bash

mus="40 120 200"

for mu in $mus ; do
  ../../main-003 --grid-size 100000 --threshold 15000 --sim 0 --topo-filename topo-$mu --mu $mu
done

rm output.csv
