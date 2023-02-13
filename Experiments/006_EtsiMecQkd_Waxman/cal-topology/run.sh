#!/bin/bash

MAIN=../main-006

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

nodes="50 100"
alpha="0.2 0.4 0.6"
beta="0.2 0.4 0.6"

for n in $nodes ; do
for a in $alpha ; do
for b in $beta ; do

  output=out-$n-$a-$b
  cmd="$MAIN \
    --output data/$output.csv \
    --num-threads $NUM_THREADS \
    --seed-start $SEED_START \
    --seed-end $SEED_END  \
    --nodes $n \
    --alpha $a \
    --beta $b \
    --max-distance 100 \
    --max-capacity 30 \
    --app-spec ../applications.dat \
    --applications 0 \
    --edge-nodes 10 \
    --edge-processing 2 \
    --algo random \
    "

  if [ "$DRY" != "" ] ; then
    echo $cmd
  else
    GLOG_v=$VERBOSE $cmd
  fi

done
done
done
