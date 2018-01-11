set term png
set output 'derive.png'
#set yrange [ -0.5 : 1.5 ]
#set xrange [ 0 : 1 ]
plot 'fort.30' title "derivé par Berstein",'fort.31' title "derive par Schema Centré"
set term x11
