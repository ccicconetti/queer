#!/bin/bash

metrics="netrate fidelity visits jain jitter"
variables="p k d q"

for m in $metrics ; do
  for v in $variables ; do
    plot.sh -k "top left" ../post/$m-*-$v.dat
    read -n 1 -p "press any key to continue..."
    killall gnuplot_x11
  done
done

