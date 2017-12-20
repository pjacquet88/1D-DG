set term png
set output 'base_l.png'
#set yrange [ -0.5 : 1.5 ]
#set xrange [ 0 : 1 ]
plot 'fort.31','fort.32','fort.33','fort.34','fort.35'
set term x11
