#!/bin/bash

if [ "$TARGET_DIR" == "" ] ; then
  echo "you must specify the directory where to save the artifacts in TARGET_DIR"
  exit 1
fi

if [ ! -d data ] ; then
  echo "data directory does not exist"
  exit 1
fi

pushd .. > /dev/null
experiment=$(basename $PWD)
popd > /dev/null
sub=$(basename $PWD)

echo $experiment-$sub.tgz

tar zcf $TARGET_DIR/$experiment-$sub.tgz data
