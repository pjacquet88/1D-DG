# animate the forward problem
plot 'Files/P'.i.'.dat' title 'P'.i

# animate the backward problem
#plot 'Files/B'.i.'.dat' title 'B'.i

# animate both
#plot 'Files/P'.i.'.dat' title 'P'.i,'Files/B'.i.'.dat' title 'B'.i

i=i+a
if (i < n) reread
