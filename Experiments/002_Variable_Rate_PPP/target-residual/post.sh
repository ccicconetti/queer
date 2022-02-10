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

residuals="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
mus="50 100"
thresholds="15000 20000"

columns=(26 27 28 29 30 31 32 33 34 35)
names=("capacity" "residual" "num-apps" "visits" "grossrate" "netrate" "pathsize" "fidelity" "jain" "jitter" )

for m in $mus ; do
  for t in $thresholds ; do
    for i in ${!columns[@]}; do
      outmangle=${names[$i]}-$m-$t
      echo "$outmangle"
      outfile=post/$outmangle.dat
      rm -f $outfile 2> /dev/null
      for r in $residuals ; do
        inmangle=out-$r-$m-$t
        datafile=data/$inmangle.csv
        value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
        echo $r $value >> $outfile
      done
    done
  done
done
