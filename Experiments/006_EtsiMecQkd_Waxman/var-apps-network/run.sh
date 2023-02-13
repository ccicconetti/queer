#!/bin/bash

MAIN=../main-006

if [ "$SEED_START" == "" ] ; then
  SEED_START=0
fi

if [ "$SEED_END" == "" ] ; then
  SEED_END=1000
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

apps="100 200 300 400 500"
algos="random spf bestfit random-blind spf-blind bestfit-blind"

for a in $apps ; do
for x in $algos ; do

  output=out-$a-$x
  cmd="$MAIN \
    --output data/$output.csv \
    --num-threads $NUM_THREADS \
    --seed-start $SEED_START \
    --seed-end $SEED_END  \
    --nodes 50 \
    --alpha 0.4 \
    --beta 0.4 \
    --max-distance 100 \
    --max-capacity 30 \
    --app-spec ../applications.dat \
    --applications $a \
    --edge-nodes 10 \
    --edge-processing 5 \
    --algo $x \
    "

  if [ "$DRY" != "" ] ; then
    echo $cmd
  else
    GLOG_v=$VERBOSE $cmd
  fi

done
done
