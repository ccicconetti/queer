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

rates="0.5 1 5 10 50"
fidelities="0.6 0.7 0.8 0.9 0.925"

for r in $rates ; do
  for f in $fidelities ; do
    output=data/out-$r-$f.csv
    cmd="$MAIN \
      --output $output \
      --num-threads $CONCURRENCY \
      --seed-start 0 \
      --seed-end 1000 \
      --link-min-epr 1 \
      --link-max-epr 400 \
      --graphml-file garr \
      --sim-duration 100 \
      --warmup-duration 10 \
      --arrival-rate 100 \
      --flow-duration 10 \
      --net-epr-rates $r \
      --q 0.5 \
      --fidelity-init 0.95 \
      --fidelity-threshold $f
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
