set term png
set output 'derive.png'
#set yrange [ -0.5 : 1.5 ]
#set xrange [ 0 : 1 ]
plot 'fort.51','fort.52'
set term x11
