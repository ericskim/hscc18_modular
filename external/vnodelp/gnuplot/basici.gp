# file basici.gp
set terminal postscript eps enh color solid "Courier" 28

set xlabel "t"
set ylabel "y_1"

set output 'lorenzi1.eps'
plot [0:6.3][-22:25]\
'lorenzi.out' u 1:($2+$3) \
	title 'lower bound' w l lt 1 lw 2,\
'lorenzi.out' u 1:($2-$3)\
	title 'upper bound' w l lt 3 lw 2

set output 'lorenzi2.eps'
set format y "%g"
set xtics 5.7,0.1,6.3
set ylabel "y_1" 0
plot [5.7:6.25][-22:0.4]\
'lorenzi.out' u 1:2:3 \
	title 'bounds'   w errorbars lw 2,\
'lorenzi.out' u 1:2   \
	title 'midpoint' w lines lt 3 lw 2 

set output 'lorenzi_excess.eps'
set ylabel "global excess" 2
set logscale y
set xtics 1
set format y "10^{%L}"
plot [0:6.3] 'lorenzi.out' u 1:4 \
	notitle w l lt 1 lw 2






