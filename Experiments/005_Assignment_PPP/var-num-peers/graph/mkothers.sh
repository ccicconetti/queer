#!/bin/bash

if [ "$METRIC" == "" ] ; then
  echo "environment variable METRIC not set"
  exit 1
fi

cp 005-vp-$METRIC-dense-10.plt 005-vp-$METRIC-dense-20.plt
cp 005-vp-$METRIC-dense-10.plt 005-vp-$METRIC-sparse-10.plt
cp 005-vp-$METRIC-dense-10.plt 005-vp-$METRIC-sparse-20.plt

sed -i -e "s/-10.dat/-20.dat/g" 005-vp-$METRIC-dense-20.plt

sed -i -e "s/dense-/sparse-/g" 005-vp-$METRIC-sparse-10.plt

sed -i -e "s/dense-/sparse-/g" 005-vp-$METRIC-sparse-20.plt
sed -i -e "s/-10.dat/-20.dat/g" 005-vp-$METRIC-sparse-20.plt
