#!/bin/bash

../../main-003 --graphml ../garr --topo-filename garr

grep -v 0,0,0 garr-0-vertices.dat > tmp.$$
mv tmp.$$ garr-0-vertices.dat

./topofilter.py < garr-0-edges.dat > tmp.$$
mv tmp.$$ garr-0-edges.dat
