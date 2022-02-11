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

ks="2 4 6 8 10 12 14 16 18 20"
quantums="1 10 20 30 40 60 80 100 200 400"
peers="2 4 6 8 10 12 14 16 18 20"
distances="1 2 3 4 5 6 7 8 9 10"

mus="50 100"
thresholds="15000 20000"

columns=(26 27 29 30 31 32 33 34 35)
names=("capacity" "residual" "visits" "grossrate" "netrate" "pathsize" "fidelity" "jain" "jitter" )

for m in $mus ; do
  for t in $thresholds ; do
    for i in ${!columns[@]}; do
      outmangle=${names[$i]}-$m-$t-k
      echo "$outmangle"
      outfile=post/$outmangle.dat
      rm -f $outfile 2> /dev/null
      d=5
      p=10
      q=10
      for k in $ks ; do
        inmangle=out-$m-$t-$d-$p-$q-$k
        datafile=data/$inmangle.csv
        value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
        echo $k $value >> $outfile
      done
    done
    for i in ${!columns[@]}; do
      outmangle=${names[$i]}-$m-$t-q
      echo "$outmangle"
      outfile=post/$outmangle.dat
      rm -f $outfile 2> /dev/null
      d=5
      p=10
      k=4
      for q in $quantums ; do
        inmangle=out-$m-$t-$d-$p-$q-$k
        datafile=data/$inmangle.csv
        value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
        echo $q $value >> $outfile
      done
    done
    for i in ${!columns[@]}; do
      outmangle=${names[$i]}-$m-$t-p
      echo "$outmangle"
      outfile=post/$outmangle.dat
      rm -f $outfile 2> /dev/null
      d=5
      q=10
      k=4
      for p in $peers ; do
        inmangle=out-$m-$t-$d-$p-$q-$k
        datafile=data/$inmangle.csv
        value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
        echo $p $value >> $outfile
      done
    done
    for i in ${!columns[@]}; do
      outmangle=${names[$i]}-$m-$t-d
      echo "$outmangle"
      outfile=post/$outmangle.dat
      rm -f $outfile 2> /dev/null
      p=10
      q=10
      k=4
      for d in $distances ; do
        inmangle=out-$m-$t-$d-$p-$q-$k
        datafile=data/$inmangle.csv
        value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
        echo $d $value >> $outfile
      done
    done
  done
done
