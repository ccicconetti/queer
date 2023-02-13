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

nodes="50 100"
capacities="15 30 100"

for n in $nodes ; do
for c in $capacities ; do

  output=out-$n-$c
  cmd="$MAIN \
    --output data/$output.csv \
    --num-threads $NUM_THREADS \
    --seed-start $SEED_START \
    --seed-end $SEED_END  \
    --nodes $n \
    --alpha 0.4 \
    --beta 0.4 \
    --max-distance 100 \
    --max-capacity $c \
    --app-spec ../applications.dat \
    --applications 1000 \
    --edge-nodes 10 \
    --edge-processing 99999 \
    --algo spf \
    "

  if [ "$DRY" != "" ] ; then
    echo $cmd
  else
    GLOG_v=$VERBOSE $cmd
  fi

done
done
