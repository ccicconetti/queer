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

columns=( 28 29 )
names=( "capacity" "residual" )

priorities="4 2 1"
perclass=("visits" "grossrate" "netrate" "pathsize" "fidelity" "jain" "jitter")
cnt=39
for value in "${perclass[@]}" ; do
  for p in $priorities ; do
    columns+=( $cnt )
    names+=( "$value-$p" )
    cnt=$(( cnt+1 ))
  done
done

if [ "$VERBOSE" != "" ] ; then
  for i in ${!columns[@]}; do
    echo "$i ${columns[$i]} ${names[$i]}"
  done
fi

topos="dense sparse"
algos="random bestfit drr"
numapps="10 20 50 100 200 500 1000"

for t in $topos ; do
for a in $algos ; do

  for i in ${!columns[@]}; do
    outmangle=${names[$i]}-$t-$a
    echo "$outmangle"
    outfile=post/$outmangle.dat
    rm -f $outfile 2> /dev/null
    for n in $numapps ; do
      datafile=data/out-$t-$a-$n.csv
      value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
      echo $n $value >> $outfile
    done
  done

done
done
