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

columns=( 28 29 30 31 32 33 34 35 36 )
names=( "capacity" "residual" "fairness-wpf" "avg-users-per-dc" "stddev-users-per-dc" "spread-users-per-dc" "avg-net-rate-per-dc" "stddev-net-rate-per-dc" "spread-net-rate-per-dc")

priorities="1"
perclass=("visits" "grossrate" "netrate" "pathsize" "fidelity" "jain" "jitter")
cnt=40
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
algos="random shortest-path load-balancing"
numapps="10 20"
numpeers="1 2 3 4 5"

for t in $topos ; do
for a in $algos ; do
for n in $numapps ; do

  for i in ${!columns[@]}; do
    outmangle=${names[$i]}-$t-$a-$n
    echo "$outmangle"
    outfile=post/$outmangle.dat
    rm -f $outfile 2> /dev/null
    for p in $numpeers ; do
      datafile=data/out-$t-$a-$n-$p.csv
      value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
      echo $p $value >> $outfile
    done
  done

done
done
done
