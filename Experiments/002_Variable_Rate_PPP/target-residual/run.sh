#!/bin/bash

MAIN=../main-002

if [[ "$DRY" == "" && ! -d "data" ]] ; then
  mkdir data 2> /dev/null
fi

if [ ! -x $MAIN ] ; then
  echo "missing executable: $MAIN"
  exit 1
fi

residuals="0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9"
mus="50 100"
thresholds="15000 20000"

for r in $residuals ; do
  for m in $mus ; do
    for t in $thresholds ; do
      output=out-$r-$m-$t

      cmd="$MAIN \
        --num-threads 100 \
        --output data/$output.csv \
        --seed-start 0 \
        --seed-end 10000 \
        --mu $m \
        --link-min-epr 1 \
        --link-max-epr 400 \
        --num-apps 10 \
        --num-peers-min 1 \
        --num-peers-max 10 \
        --distance-min 1 \
        --distance-max 5 \
        --grid-size 60000 \
        --threshold $t \
        --link-probability 0.5 \
        --q 0.5 \
        --quantum 10 \
        --k 4 \
        --fidelity-init 0.95 \
        --fidelity-threshold 0 \
        --target-residual $r"

      if [ "$DRY" != "" ] ; then
        echo $cmd
      else
        GLOG_v=$VERBOSE $cmd
      fi
    done
  done
done
