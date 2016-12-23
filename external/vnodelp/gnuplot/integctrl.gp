# file integctrl.gp
set terminal postscript eps enh color solid "Courier" 28

set ylabel "y_1"
set xlabel "t"

set output 'lorenz.eps'
plot [][-20:22]\
 'lorenz.tight' u 1:2 title 'lower bound' w l lw 2,\
 'lorenz.tight' u 1:3 title 'upper bound' w l lt 3 lw 2

set output 'lorenz2.eps'
set xtics 0, 0.1, 0.4
plot [0:0.4]\
 'lorenz.tight'   u 1:2 title 'tight'      w l lt 3 lw 2,\
 'lorenz.tight'   u 1:3 notitle            w l lt 3 lw 2,\
 'lorenz.apriori' u 1:2 title   'a priori' w l lt 1 lw 3

set output 'lorenz_err.eps'
set xtics 0,2,20
set ylabel "global excess" 2
set logscale y
set format y "10^{%L}"
plot 'lorenz.tight' u 1:4 notitle w l lw 2

set output 'lorenz_step.eps'
set nologscale
set ylabel "stepsize"
plot 'lorenz.step' u 1:2 notitle w l lw 2

