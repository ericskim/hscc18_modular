# file odeparam.gp
set terminal postscript eps enh color solid "Courier" 28

set xlabel "y_1"
set ylabel "y_2"
set zlabel "y_3"

set xtics -20,10,20
set ytics -25,10,25
set ztics  0,10,45

set output 'odeparam.eps'
splot 'odeparam1.out' u 1:2:3 title '{/Symbol b} = 8/3'\
		 w l lt 1 lw 2,\
      'odeparam2.out' u 1:2:3 title '{/Symbol b} = 5'\
		 w l lt 3 lw 2




