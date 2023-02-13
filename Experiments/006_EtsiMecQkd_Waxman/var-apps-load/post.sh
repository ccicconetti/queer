#!/bin/bash

percentile_script=$(which percentile.py)
if [ "$percentile_script" == "" ] ; then
  ../../../Scripts/download-percentile.sh
  percentile_script=./percentile.py
fi

if [ ! -d "data" ] ; then
  echo "no 'data' directory"
fi

if [ ! -d "post" ] ; then
  mkdir post 2> /dev/null
fi

columns=( 19 20 21 22 23 24 25 26 27 28 )
names=( "capacity-tot" "processing-tot" "capacity-res" "processing-res" "num-apps-allocated" "path-length-avg" "net-rate-tot" "edge-node-util-stddev" "edge-node-util-jain" "edge-node-util-spread" )

if [ "$VERBOSE" != "" ] ; then
  for i in ${!columns[@]}; do
    echo "$i ${columns[$i]} ${names[$i]}"
  done
fi

apps="100 200 300 400 500 600 700"
algos="random spf bestfit random-blind spf-blind bestfit-blind"

for x in $algos ; do

  for i in ${!columns[@]}; do
    outmangle=${names[$i]}-$x
    echo "$outmangle"
    outfile=post/$outmangle.dat
    rm -f $outfile 2> /dev/null
    for a in $apps ; do
      datafile=data/out-$a-$x.csv
      value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
      echo $a $value >> $outfile
    done
  done

done
