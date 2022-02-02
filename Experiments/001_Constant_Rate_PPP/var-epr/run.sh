#!/bin/bash

MAIN=../main-001

if [[ "$DRY" == "" && ! -d "data" ]] ; then
  mkdir data 2> /dev/null
fi

if [ ! -x $MAIN ] ; then
  echo "missing executable: $MAIN"
  exit 1
fi

flows="100 1000 10000"
eprs="49 99 149 199 249 299 349 399 449 499"

for f in $flows ; do
  for e in $eprs ; do

    cmd="$MAIN \
      --num-threads 100 \
      --output data/out-$f-$e.csv \
      --seed-start 0 \
      --seed-end 10000 \
      --mu 100 \
      --num-flows $f \
      --link-min-epr 1 \
      --link-max-epr $e \
      --rate-min 1 \
      --rate-max 10"

    if [ "$DRY" != "" ] ; then
      echo $cmd
    else
      GLOG_v=$VERBOSE $cmd
    fi
  done
done
