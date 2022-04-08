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

numnodes="50 75 100 125 150 175 200"

for n in $numnodes ; do
  output=data/out-$n.csv
  cmd="$MAIN \
    --output $output \
    --num-threads $CONCURRENCY \
    --mu $n \
    --seed-start 0 \
    --seed-end 64 \
    --link-min-epr 1 \
    --link-max-epr 1 \
    --sim-duration 0 \
    --warmup-duration 0 \
    --arrival-rate 1 \
    --flow-duration 10 \
    --net-epr-rates 1 \
    --grid-size 100000 \
    --threshold 15000 \
    --link-probability 1 \
    --q 0.5 \
    --fidelity-init 0.95 \
    --fidelity-threshold 1 \
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
