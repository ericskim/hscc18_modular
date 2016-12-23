# file orbit.gp
set terminal postscript eps enh color solid "Courier" 28

set xlabel "y_1"
set ylabel "y_2"


set output 'orbit_sol.eps'
plot  'orbit_sol.out' u 1:2 notitle w l lt 1 lw 2

set xlabel "t"
set ylabel "stepsize"
set output 'orbit_step.eps'
plot  'orbit_step.out' u 1:2 notitle w l lt 2 lw 2

