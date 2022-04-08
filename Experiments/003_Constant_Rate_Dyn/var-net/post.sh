#!/bin/bash

if [ ! -d "data" ] ; then
  echo "no 'data' directory"
fi

if [ ! -d "post" ] ; then
  mkdir post 2> /dev/null
fi

columns=(24 25 26 30 28 29 31 32)
names=("capacity" "residual" "num-active-flows" "admission-rate" "gross-rate" "net-rate" "path-size" "fidelity")

maxlinkrates="20 40 60 80 100"
numnodes="40 80 120 160 200"
rates="1 10"
fidelities="0.7 0.9"

for i in ${!columns[@]}; do
  for r in $rates ; do
    for f in $fidelities ; do
      for m in $maxlinkrates ; do
        for n in $numnodes ; do
          outmangle=${names[$i]}-$r-$f
          datafile=data/out-$m-$n-$r-$f.csv
          echo "$outmangle"

          outfile=post/$outmangle.dat
          rm -f $outfile 2> /dev/null
          cut -d , -f 2,7,${columns[$i]} < $datafile  | tr ',' ' ' >> $outfile
          echo >> $outfile
        done
      done
    done
  done
done
