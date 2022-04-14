#!/bin/bash

if [ -d data ] ; then
  echo "data directory already present: cowardly refusing to proceed"
  exit 1
fi

pushd .. > /dev/null
experiment=$(basename $PWD)
popd > /dev/null
sub=$(basename $PWD)

mkdir data 2> /dev/null
wget http://turig.iit.cnr.it/~claudio/public/quantum-routing/$experiment-$sub.tgz -O- | tar zx
