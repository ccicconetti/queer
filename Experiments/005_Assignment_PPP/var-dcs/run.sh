#!/bin/bash

MAIN=../main-005

if [ "$SEED_START" == "" ] ; then
  SEED_START=0
fi

if [ "$SEED_END" == "" ] ; then
  SEED_END=5000
fi

if [ "$NUM_THREADS" == "" ] ; then
  NUM_THREADS=0
fi

if [[ "$DRY" == "" && ! -d "data" ]] ; then
  mkdir data 2> /dev/null
fi

if [ ! -x $MAIN ] ; then
  echo "missing executable: $MAIN"
  exit 1
fi

topos="dense sparse"
algos="random shortest-path load-balancing"
fracs="0.1 0.2 0.3 0.4 0.5"

for t in $topos ; do
for a in $algos ; do
for f in $fracs ; do

  if [ $t == "dense" ] ; then
    mu=100
    threshold=20000
  elif [ $t == "sparse" ] ; then
    mu=50
    threshold=15000
  else
    echo "invalid topology type: $t"
    exit 1
  fi

  output=out-$t-$a-$f
  cmd="$MAIN \
    --num-threads $NUM_THREADS \
    --seed-start $SEED_START \
    --seed-end $SEED_END  \
    --mu $mu \
    --link-min-epr 1 \
    --link-max-epr 400 \
    --num-apps 15 \
    --num-peers 3 \
    --peer-assignment-algo $a \
    --grid-size 60000 \
    --threshold $threshold \
    --link-probability 0.5 \
    --algorithm drr \
    --q 0.5 \
    --quantum 10 \
    --k 4 \
    --fidelity-init 0.95 \
    --frac-end-users 0.2 \
    --frac-data-centers $f \
    --priorities 1 \
    --fidelity-thresholds 0.5 \
    --target-residual -1 \
    --output data/$output.csv \
    "

  if [ "$DRY" != "" ] ; then
    echo $cmd
  else
    GLOG_v=$VERBOSE $cmd
  fi

done
done
done
