#!/bin/bash

MAIN=../main-003

if [[ "$DRY" == "" && ! -d "data" ]] ; then
  mkdir data 2> /dev/null
fi

if [[ "$CONCURRENCY" == "" ]] ; then
  if [ -r /proc/cpuinfo ] ; then
    CONCURRENCY=$(cat /proc/cpuinfo  | grep "core id" | wc -l)
  else
    CONCURRENCY=1
  fi
fi

if [ ! -x $MAIN ] ; then
  echo "missing executable: $MAIN"
  exit 1
fi

srcdstpolicies="uniform nodecapacities"
maxlinkrates="100 150 200 250 300 350 400 450 500"
numnodes="40 60 80 100 120 140 160 180 200"
rates="1 10"
fidelities="0.7 0.9"

for p in $srcdstpolicies ; do
  for m in $maxlinkrates ; do
    for n in $numnodes ; do
      for r in $rates ; do
        for f in $fidelities ; do
          output=data/out-$p-$m-$n-$r-$f.csv
          cmd="$MAIN \
            --output $output \
            --num-threads $CONCURRENCY \
            --seed-start 0 \
            --seed-end 100 \
            --mu $n \
            --link-min-epr 1 \
            --link-max-epr $m \
            --src-dst-policy $p \
            --grid-size 100000 \
            --threshold 15000 \
            --sim-duration 100 \
            --warmup-duration 10 \
            --arrival-rate 100 \
            --flow-duration 10 \
            --net-epr-rates $r \
            --q 0.5 \
            --fidelity-init 0.95 \
            --fidelity-threshold $f \
            "

          if [ "$EXPLAIN" != "" ] ; then
            $cmd --explain
            exit 1

          elif [ "$DRY" != "" ] ; then
            echo $cmd

          else
            now=$(date)
            echo -n "$now $output."
            if [ -r $output ] ; then
              echo ".skipped"
            else
              GLOG_v=$VERBOSE $cmd
              echo ".done"
            fi
          fi
        done
      done
    done
  done
done
