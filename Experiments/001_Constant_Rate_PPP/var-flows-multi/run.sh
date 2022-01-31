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
mus="50 100"
eprs="constant uniform"

for f in $flows ; do
  for m in $mus ; do
    for e in $eprs ; do

      if [ $e == "constant" ] ; then
        link_epr="--link-min-epr 1 --link-max-epr 99"
      elif [ $e == "uniform" ] ; then
        link_epr="--link-min-epr 50 --link-max-epr 50"
      else
        echo "invalid link EPR: $e"
        exit 1
      fi

      cmd="$MAIN \
        --num-threads 100 \
        --output data/out-$f-$m-$e.csv \
        --seed-start 0 \
        --seed-end 999 \
        --mu $m \
        --num-flows $f \
        $link_epr \
        --rate-min 1 \
        --rate-max 10"

      if [ "$DRY" != "" ] ; then
        echo $cmd
      else
        GLOG_v=$VERBOSE $cmd
      fi
    done
  done
done
