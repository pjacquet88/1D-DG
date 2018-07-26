# animate the forward problem
#plot 'Files/P'.i.'.dat' title 'P'.i,'Files/U'.i.'.dat' title 'U'.i

# animate the backward problem
#plot 'Files/B'.i.'.dat' title 'B'.i

# animate both
#plot 'Files/P'.i.'.dat' title 'P'.i,'Files/B'.i.'.dat' title 'B'.i


# test
# plot 'Files/FU'.i.'.dat' title 'FU'.i,'Files/U'.i.'.dat' title 'U'.i,'Files/FP'.i.'.dat' title 'FP'.i,'Files/P'.i.'.dat' title 'P'.i
# plot 'Files/FU'.i.'.dat' title 'FU'.i,'Files/FP'.i.'.dat' title 'FP'.i
plot 'Files/FP'.i.'.dat' title 'FP'.i ,'Files/P'.i.'.dat' title 'P'.i

i=i+a
if (i < n) reread
