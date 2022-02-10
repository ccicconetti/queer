#!/bin/bash

MAIN=../main-002

if [[ "$DRY" == "" && ! -d "data" ]] ; then
  mkdir data 2> /dev/null
fi

if [ ! -x $MAIN ] ; then
  echo "missing executable: $MAIN"
  exit 1
fi

mus="50 100"
linkthresholds="15000 20000"
linkprobabilities="0.5 0.75 1"

for m in $mus ; do
  for t in $linkthresholds ; do
    for p in $linkprobabilities ; do
      output="out-$m-$t-$p"

      cmd="$MAIN \
        --num-threads 100 \
        --output data/$output.csv \
        --seed-start 0 \
        --seed-end 1000 \
        --mu $m \
        --num-apps 0 \
        --link-probability $p \
        --threshold $t"

      if [ "$DRY" != "" ] ; then
        echo $cmd
      else
        GLOG_v=$VERBOSE $cmd
      fi
    done
  done
done
