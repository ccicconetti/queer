#!/bin/bash

MAIN=../main-007

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

arrival_rates="0.0001 0.0002 0.0005 0.001 0.002 0.005"
algos="policy-014-k-1 policy-014-k-3 policy-015 policy-015-reuse"

for a in $arrival_rates ; do
for x in $algos ; do

  output=out-$a-$x
  cmd="$MAIN \
    --output data/$output.csv \
    --num-threads $NUM_THREADS \
    --seed-start $SEED_START \
    --seed-end $SEED_END  \
    --duration $((86400*14)) \
    --warmup $((86400*2)) \
    --nodes 50 \
    --alpha 0.4 \
    --beta 0.4 \
    --max-distance 100 \
    --max-capacity 30 \
    --app-spec ../applications.dat \
    --arrival-rate $a \
    --edge-nodes 10 \
    --edge-processing U(3,7) \
    --algo $x \
    "

  if [ "$DRY" != "" ] ; then
    echo $cmd
  else
    GLOG_v=$VERBOSE $cmd
  fi

done
done