set term png
set output 'base_b.png'
#set yrange [ -0.5 : 1.5 ]
#set xrange [ 0 : 1 ]
plot for [i=101:130] 'fort.'.i
set term x11
