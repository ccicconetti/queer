set term x11 persist
set grid
set key top left
set xlabel "Average number of peers/host"
set ylabel "Jain's fairness index"
plot '../post/jain-100-15000-p.dat' u 1:2 w lp pt 6 lt 1 title "{/Symbol m} = 100, {/Symbol t} = 15 km",'../post/jain-100-20000-p.dat' u 1:2 w lp pt 7 lt 2 title "{/Symbol m} = 100, {/Symbol t} = 20 km",'../post/jain-50-15000-p.dat' u 1:2 w lp pt 8 lt 3 title "{/Symbol m} = 50, {/Symbol t} = 15 km",'../post/jain-50-20000-p.dat' u 1:2 w lp pt 9 lt 4 title "{/Symbol m} = 50, {/Symbol t} = 20 km"
