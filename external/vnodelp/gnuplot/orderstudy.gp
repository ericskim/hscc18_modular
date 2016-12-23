# file orderstudy.gp
set terminal postscript eps enh color solid "Courier" 28

set ylabel "CPU time"
set xlabel "order"

set output 'order.eps'
plot  [10:35]\
      'order1e-07.out' u 1:2 title '10^{-7}'  w l lt 1 lw 2,\
      'order1e-08.out' u 1:2 title '10^{-8}'  w l lt 2 lw 2,\
      'order1e-09.out' u 1:2 title '10^{-9}'  w l lt 3 lw 2,\
      'order1e-10.out' u 1:2 title '10^{-10}' w l lt 4 lw 2,\
      'order1e-11.out' u 1:2 title '10^{-11}' w l lt 5 lw 2,\
      'order1e-12.out' u 1:2 title '10^{-12}' w l lt 9 lw 2,\
      'order1e-13.out' u 1:2 title '10^{-13}' w l lt 7 lw 2
     

set logscale y
set format y "10^{%L}"
set output 'timeorder.eps'
plot 'order1e-07.out' u 1:2 title '10^{-7}'  w l lt 1 lw 2,\
     'order1e-08.out' u 1:2 title '10^{-8}'  w l lt 2 lw 2,\
     'order1e-09.out' u 1:2 title '10^{-9}'  w l lt 3 lw 2,\
     'order1e-10.out' u 1:2 title '10^{-10}' w l lt 4 lw 2,\
     'order1e-11.out' u 1:2 title '10^{-11}' w l lt 5 lw 2,\
     'order1e-12.out' u 1:2 title '10^{-12}' w l lt 9 lw 2,\
     'order1e-13.out' u 1:2 title '10^{-13}' w l lt 7 lw 2




