#!/bin/bash

MAIN=../main-004

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
algos="random bestfit drr"
numapps="10 20 50 100 200 500 1000"

for t in $topos ; do
for a in $algos ; do
for n in $numapps ; do

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

  output=out-$t-$a-$n
  cmd="$MAIN \
    --num-threads $NUM_THREADS \
    --seed-start $SEED_START \
    --seed-end $SEED_END  \
    --mu $mu \
    --link-min-epr 1 \
    --link-max-epr 400 \
    --num-apps $n \
    --num-peers-min 2 \
    --num-peers-max 4 \
    --distance-min 2 \
    --distance-max 7 \
    --grid-size 60000 \
    --threshold $threshold \
    --link-probability 0.5 \
    --algorithm $a \
    --q 0.5 \
    --quantum 10 \
    --k 4 \
    --fidelity-init 0.95 \
    --priorities 1 \
    --fidelity-thresholds 0.7,0.8,0.9 \
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
