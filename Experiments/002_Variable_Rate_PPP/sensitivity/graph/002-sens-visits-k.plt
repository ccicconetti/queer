#set term x11 persist
set grid
set key top left
set xlabel "Maximum number of shortest paths per peer [k]"
set ylabel "Average number of visits"
plot '../post/visits-100-15000-k.dat' u 1:2 w lp pt 6 lt 1 title "{/Symbol m} = 100, {/Symbol t} = 15 km",'../post/visits-100-20000-k.dat' u 1:2 w lp pt 7 lt 2 title "{/Symbol m} = 100, {/Symbol t} = 20 km",'../post/visits-50-15000-k.dat' u 1:2 w lp pt 8 lt 3 title "{/Symbol m} = 50, {/Symbol t} = 15 km",'../post/visits-50-20000-k.dat' u 1:2 w lp pt 9 lt 4 title "{/Symbol m} = 50, {/Symbol t} = 20 km"
