# file work.gp
set terminal postscript eps enh color solid "Courier" 28

# model of the work	
f(x) = a + b*x

fit f(x) "work.out" using (log($1)):(log($2)) via a,b

set xlabel "number equations"
set ylabel "time per step (seconds)"
set xrange [40:200]	

set output 'work.eps'	
plot 'work.out' using 1:2 notitle with lines

set logscale
set output 'worklog.eps'
plot 'work.out' using 1:2 notitle with lines
