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

columns=(31 26 27 25 29 30 32 33)
names=("admission-rate" "residual" "num-active-flows" "capacity" "gross-rate" "net-rate" "path-size" "fidelity")

qvalues="0.5 0.6 0.7 0.8 0.9 1"
fidelities="0.95 0.97 0.99 0.999"

for f in $fidelities ; do
  for i in ${!columns[@]}; do
    outmangle=${names[$i]}-$f
    outfile=post/$outmangle.dat
    rm -f $outfile 2> /dev/null
    echo "$outmangle"
    for q in $qvalues ; do
      datafile=data/out-$q-$f.csv
      value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
      echo "$q $value" >> $outfile
    done
  done
done
