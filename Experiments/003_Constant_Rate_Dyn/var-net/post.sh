#!/bin/bash

if [ ! -d "data" ] ; then
  echo "no 'data' directory"
fi

percentile_script=$(which percentile.py)
if [ "$percentile_script" == "" ] ; then
  if [ ! -x percentile.py ] ; then
    curl -opercentile.py https://raw.githubusercontent.com/ccicconetti/serverlessonedge/master/scripts/percentile.py >& /dev/null
    if [ $? -ne 0 ] ; then
      echo "error downloading the percentile.py script"
    fi
    chmod 755 percentile.py
  fi
  percentile_script=./percentile.py
fi


if [ ! -d "post" ] ; then
  mkdir post 2> /dev/null
fi

columns=(30 25 26 24 28 29 31 32)
names=("admission-rate" "residual" "num-active-flows" "capacity" "gross-rate" "net-rate" "path-size" "fidelity")

maxlinkrates="20 30 40 50 60 70 80 90 100"
numnodes="40 60 80 100 120 140 160 180 200"
rates="1 10"
fidelities="0.7 0.9"

for i in ${!columns[@]}; do
  for r in $rates ; do
    for f in $fidelities ; do
      outmangle=${names[$i]}-$r-$f
      outfile=post/$outmangle.dat
      rm -f $outfile 2> /dev/null
      echo "$outmangle"
      for m in $maxlinkrates ; do
        for n in $numnodes ; do
          datafile=data/out-$m-$n-$r-$f.csv
          value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
          echo "$m $n $value" >> $outfile
        done
        echo >> $outfile
      done
    done
  done
done
