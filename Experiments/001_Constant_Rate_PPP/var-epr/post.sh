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

flows="100 1000 10000"
eprs="49 99 149 199 249 299 349 399 449 499"
columns=(14 15 20 21 22 23 24 25 27 28)
names=("num-nodes" "num-edges" "capacity" "residual" "dijkstracalls" "grossrate" "netrate" "admission" "pathsize" "fidelity")

for f in $flows ; do
  for i in ${!columns[@]}; do
    outmangle=${names[$i]}-$f
    echo "$outmangle"
    outfile=post/$outmangle.dat
    rm -f $outfile 2> /dev/null
    for e in $eprs ; do
      inmangle=out-$f-$e
      datafile=data/$inmangle.csv
      value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
      echo $e $value >> $outfile
    done
  done
done
