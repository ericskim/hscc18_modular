# file vanderpol.gp
set terminal postscript eps enh color solid "Courier" 28

set xlabel 't' 
set ylabel 'stepsize'

set output 'vdp_step1.eps'
plot [0:200][0:0.18] 'vdp_step1.out' u 1:2 \
	title '{/Symbol m}=10^1' w l lt 2 lw 2,\
     'vdp_step2.out' u 1:2 \
	title '{/Symbol m}=10^2' w l lt 3 lw 2

set output 'vdp_step2.eps'
plot [0:200][0:0.002] 'vdp_step3.out' u 1:2 \
	title '{/Symbol m}=10^3' w l lt 3 lw 2,\
     'vdp_step4.out' u 1:2 \
	title '{/Symbol m}=10^4' w l lt 4 lw 2
	





