#!/bin/bash

MAIN=../main-001

if [[ "$DRY" == "" && ! -d "data" ]] ; then
  mkdir data 2> /dev/null
fi

if [ ! -x $MAIN ] ; then
  echo "missing executable: $MAIN"
  exit 1
fi

flows="10 20 50 100 200 500 1000 2000 5000 10000"

for f in $flows ; do
  cmd="$MAIN \
    --num-threads 100 \
    --output data/out-$f.csv \
    --seed-start 0 \
    --seed-end 10000 \
    --mu 50 \
    --link-min-epr 1 \
    --link-max-epr 100 \
    --num-flows $f \
    --rate-min 1 \
    --rate-max 10"

  if [ "$DRY" != "" ] ; then
    echo $cmd
  else
    GLOG_v=$VERBOSE $cmd
  fi
done
