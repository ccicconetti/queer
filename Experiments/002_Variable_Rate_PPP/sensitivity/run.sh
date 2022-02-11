#!/bin/bash

MAIN=../main-002

if [[ "$DRY" == "" && ! -d "data" ]] ; then
  mkdir data 2> /dev/null
fi

if [ ! -x $MAIN ] ; then
  echo "missing executable: $MAIN"
  exit 1
fi

ks="2 4 6 8 10 12 14 16 18 20"
quantums="1 10 20 30 40 60 80 100 200 400"
peers="2 4 6 8 10 12 14 16 18 20"
distances="1 2 3 4 5 6 7 8 9 10"

mus="50 100"
thresholds="15000 20000"

for m in $mus ; do
  for t in $thresholds ; do
    cmdbase="$MAIN \
      --num-threads 100 \
      --seed-start 0 \
      --seed-end 1000 \
      --mu $m \
      --link-min-epr 1 \
      --link-max-epr 400 \
      --num-apps 100 \
      --num-peers-min 1 \
      --distance-min 1 \
      --grid-size 60000 \
      --threshold $t \
      --link-probability 0.5 \
      --q 0.5 \
      --fidelity-init 0.95 \
      --fidelity-threshold 0 \
      --target-residual -1"

    d=5
    p=10
    q=10
    for k in $ks ; do
      output=out-$m-$t-$d-$p-$q-$k
      cmd="$cmdbase \
        --output data/$output.csv \
        --distance-max $d \
        --num-peers-max $p \
        --quantum $q \
        --k $k"
      if [ "$DRY" != "" ] ; then
        echo $cmd
      else
        GLOG_v=$VERBOSE $cmd
      fi
    done

    d=5
    p=10
    k=4
    for q in $quantums ; do
      output=out-$m-$t-$d-$p-$q-$k
      cmd="$cmdbase \
        --output data/$output.csv \
        --distance-max $d \
        --num-peers-max $p \
        --quantum $q \
        --k $k"
      if [ "$DRY" != "" ] ; then
        echo $cmd
      else
        GLOG_v=$VERBOSE $cmd
      fi
    done

    d=5
    q=10
    k=4
    for p in $peers ; do
      output=out-$m-$t-$d-$p-$q-$k
      cmd="$cmdbase \
        --output data/$output.csv \
        --distance-max $d \
        --num-peers-max $p \
        --quantum $q \
        --k $k"
      if [ "$DRY" != "" ] ; then
        echo $cmd
      else
        GLOG_v=$VERBOSE $cmd
      fi
    done

    p=10
    q=10
    k=4
    for d in $distances ; do
      output=out-$m-$t-$d-$p-$q-$k
      cmd="$cmdbase \
        --output data/$output.csv \
        --distance-max $d \
        --num-peers-max $p \
        --quantum $q \
        --k $k"
      if [ "$DRY" != "" ] ; then
        echo $cmd
      else
        GLOG_v=$VERBOSE $cmd
      fi
    done

  done
done
