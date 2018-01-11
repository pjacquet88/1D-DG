set term png
set output 'b2l.png'
#set yrange [ -0.5 : 1.5 ]
#set xrange [ 0 : 1 ]
plot 'fort.20' title "Bernstein",'fort.21' title "Lagrange"
set term x11
