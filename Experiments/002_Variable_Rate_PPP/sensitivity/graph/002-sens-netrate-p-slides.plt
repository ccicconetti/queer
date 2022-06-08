#set term x11 persist
set grid
set key top left
set xlabel "Maximum number of peers/host"
set ylabel "Net rate (Bell pairs/s)"
set yrange [0:*]
unset key
unset grid
plot '../post/netrate-100-15000-p.dat' u 1:2 w l lw 3
