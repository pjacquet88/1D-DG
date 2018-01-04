set term png
set output 'l2b.png'
#set yrange [ -0.5 : 1.5 ]
#set xrange [ 0 : 1 ]
plot 'fort.151','fort.152'
set term x11
