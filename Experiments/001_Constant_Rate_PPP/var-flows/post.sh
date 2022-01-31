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
columns=(20 21 22 23 24 25 27)
names=("capacity" "residual" "dijkstracalls" "grossrate" "netrate" "admission" "pathsize")

for i in ${!columns[@]}; do
  echo "${names[$i]}"
  outfile=post/${names[$i]}.dat
  rm -f $outfile 2> /dev/null
  for f in $flows ; do
    datafile=data/out-$f.csv
    value=$(sed -e "s/,/ /g" $datafile | percentile.py --column ${columns[$i]} --mean | cut -f 1 -d ' ')
    echo $f $value >> $outfile
  done
done
