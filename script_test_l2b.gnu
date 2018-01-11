set term png
set output 'l2b.png'
#set yrange [ -0.5 : 1.5 ]
#set xrange [ 0 : 1 ]
plot 'fort.10' title "Bernstein" ,'fort.11' title "Lagrange"
set term x11
