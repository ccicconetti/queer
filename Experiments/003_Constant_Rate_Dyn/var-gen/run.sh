#!/bin/bash

MAIN=../main-003

if [[ "$DRY" == "" && ! -d "data" ]] ; then
  mkdir data 2> /dev/null
fi

if [[ "$CONCURRENCY" == "" ]] ; then
  if [ -r /proc/cpuinfo ] ; then
    CONCURRENCY=$(cat /proc/cpuinfo  | grep "core id" | wc -l)
  else
    CONCURRENCY=1
  fi
fi

if [ ! -x $MAIN ] ; then
  echo "missing executable: $MAIN"
  exit 1
fi

qvalues="0.5 0.6 0.7 0.8 0.9 1"
fidelities="0.95 0.97 0.99 0.999"

for q in $qvalues ; do
  for f in $fidelities ; do
    output=data/out-$q-$f.csv
    cmd="$MAIN \
      --output $output \
      --num-threads $CONCURRENCY \
      --seed-start 0 \
      --seed-end 100 \
      --mu 120 \
      --link-min-epr 1 \
      --link-max-epr 300 \
      --src-dst-policy uniform \
      --grid-size 100000 \
      --threshold 15000 \
      --sim-duration 1000 \
      --warmup-duration 50 \
      --arrival-rate 100 \
      --flow-duration 10 \
      --net-epr-rates 10 \
      --q $q \
      --fidelity-init $f \
      --fidelity-threshold 0.9 \
      "

    if [ "$EXPLAIN" != "" ] ; then
      $cmd --explain
      exit 1

    elif [ "$DRY" != "" ] ; then
      echo $cmd

    else
      now=$(date)
      echo -n "$now $output."
      if [ -r $output ] ; then
        echo ".skipped"
      else
        GLOG_v=$VERBOSE $cmd
        echo ".done"
      fi
    fi
  done
done
