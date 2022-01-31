#!/bin/bash

if [ ! -d "data" ] ; then
  echo "no 'data' directory"
fi

if [ "$(which percentile.py)" == "" ] ; then
  echo "wrong environment"
  exit 1
fi

if [ ! -d "post" ] ; then
  mkdir post 2> /dev/null
fi

flows="10 20 50 100 200 500 1000 2000 5000 10000"
mus="50 100"
eprs="constant uniform"
columns=(20 21 22 23 24 25 27)
names=("capacity" "residual" "dijkstracalls" "grossrate" "netrate" "admission" "pathsize")

for m in $mus ; do
  for e in $eprs ; do
    for i in ${!columns[@]}; do
      outmangle=${names[$i]}-$m-$e
      echo "$outmangle"
      outfile=post/$outmangle.dat
      rm -f $outfile 2> /dev/null
      for f in $flows ; do
        inmangle=out-$f-$m-$e
        datafile=data/$inmangle.csv
        value=$(sed -e "s/,/ /g" $datafile | percentile.py --column ${columns[$i]} --mean | cut -f 1,3 -d ' ')
        echo $f $value >> $outfile
      done
    done
  done
done
