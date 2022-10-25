#!/bin/bash

if [ "$SERVERS" == "" ] ; then
  echo "empty environment variable SERVERS"
  exit 1
fi

if [ -d data ] ; then
  read -n 1 -p "Do you want to remove the directory 'data'? (Ctrl+C to abort)"
  rm -rf data >& /dev/null
fi

for i in $SERVERS ; do
  scp -r $i:$PWD/data data-$i
done

mkdir data 2> /dev/null
for i in $(cd data-$(echo $SERVERS | cut -f 1 -d ' ') && ls -1) ; do
  cat data-*/$i > data/$i
done
