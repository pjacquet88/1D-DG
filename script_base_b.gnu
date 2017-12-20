set term png
set output 'base_b.png'
#set yrange [ -0.5 : 1.5 ]
#set xrange [ 0 : 1 ]
plot 'fort.21','fort.22','fort.23','fort.24','fort.25'
set term x11
