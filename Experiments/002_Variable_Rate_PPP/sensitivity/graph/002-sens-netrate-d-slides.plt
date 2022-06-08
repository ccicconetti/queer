#set term x11 persist
set grid
set key top right
set xlabel "Maximum distance of peers (hops)"
set ylabel "Net rate (Bell pairs/s)"
unset grid
unset key
plot '../post/netrate-100-15000-d.dat' u 1:2 w l lw 3
