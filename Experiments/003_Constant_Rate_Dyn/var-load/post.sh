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

columns=(24 25 26 30 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48)
names=("capacity" "residual" "num-active-flows" "admission-rate" "gross-rate-1-0.7" "gross-rate-1-0.9" "gross-rate-10-0.7" "gross-rate-10-0.9" "net-rate-1-0.7" "net-rate-1-0.9" "net-rate-10-0.7" "net-rate-10-0.9" "admission-rate-1-0.7" "admission-rate-1-0.9" "admission-rate-10-0.7" "admission-rate-10-0.9" "avg-path-size-1-0.7" "avg-path-size-1-0.9" "avg-path-size-10-0.7" "avg-path-size-10-0.9")

arrivalrates="1 5 10 50 100 500 1000"
graphmls="garr"

for g in $graphmls ; do
  for i in ${!columns[@]}; do
    outmangle=${names[$i]}-$g
    echo "$outmangle"
    outfile=post/$outmangle.dat
    rm -f $outfile 2> /dev/null
    for a in $arrivalrates ; do
      inmangle=out-$g-$a
      datafile=data/$inmangle.csv
      value=$($percentile_script --delimiter , --column ${columns[$i]} --mean < $datafile | cut -f 1,3 -d ' ')
      echo $a $value >> $outfile
    done
  done
done
