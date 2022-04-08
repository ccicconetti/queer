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

arrivalrates="1 5 10 50 100 500 1000"
graphmls="garr"

for g in $graphmls ; do
  for a in $arrivalrates ; do
    output=data/out-$g-$a.csv
    cmd="$MAIN \
      --output $output \
      --num-threads $CONCURRENCY \
      --seed-start 0 \
      --seed-end 10000 \
      --link-min-epr 1 \
      --link-max-epr 400 \
      --graphml-file $g \
      --sim-duration $(echo "10000/$a" | bc) \
      --warmup-duration $(echo "1000/$a" | bc) \
      --arrival-rate $a \
      --flow-duration 10 \
      --net-epr-rates 1@10 \
      --q 0.5 \
      --fidelity-init 0.95 \
      --fidelity-threshold 0.7@0.9 \
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
