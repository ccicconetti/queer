#!/bin/bash

metrics=("netrate" "fidelity" "visits" "jain" "jitter")
metrics_labels=("Net rate (EPR-pairs/s)" "Fidelity" "Average number of visits" "Jain's fairness index" "Max-min net rate (EPR-pairs/s)")
variables=("p" "k" "d" "q")
variables_labels=("Average number of peers/host" "Maximum number of shortest paths per peer [k]" "Maximum distance of peers" "Scheduler quantum (EPR-pairs/s)")

for m in ${!metrics[@]} ; do
  for v in ${!variables[@]} ; do
    metric=${metrics[$m]}
    variable=${variables[$v]}
    plotfile=$metric-$variable.plt

    if [ -r $plotfile ] ; then
      echo "skipping: $plotfile"
    else
      echo "writing to: $plotfile"
    fi

    metric_label=${metrics_labels[$m]}
    variable_label=${variables_labels[$v]}
    cat << EOF > $plotfile
set term x11 persist
set grid
set key top left
set xlabel "$variable_label"
set ylabel "$metric_label"
plot \
'../post/$metric-100-15000-$variable.dat' u 1:2 w lp pt 6 lt 1 title "{/Symbol m} = 100, {/Symbol t} = 15 km",\
'../post/$metric-100-20000-$variable.dat' u 1:2 w lp pt 7 lt 2 title "{/Symbol m} = 100, {/Symbol t} = 20 km",\
'../post/$metric-50-15000-$variable.dat' u 1:2 w lp pt 8 lt 3 title "{/Symbol m} = 50, {/Symbol t} = 15 km",\
'../post/$metric-50-20000-$variable.dat' u 1:2 w lp pt 9 lt 4 title "{/Symbol m} = 50, {/Symbol t} = 20 km"
EOF
  done
done

