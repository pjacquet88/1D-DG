# animate the forward problem
#plot 'Files/P'.i.'.dat' title 'P'.i#,'Files/QP'.i.'.dat' title 'QP'.i
#plot 'Files/FP'.i.'.dat' title 'FP'.i
#plot 'Files/QP'.i.'.dat' title 'QP'.i with lines

# animate the backward problem
#plot 'Files/B'.i.'.dat' title 'B'.i

# animate both
#plot 'Files/P'.i.'.dat' title 'P'.i,'Files/B'.i.'.dat' title 'B'.i


# test
# plot 'Files/FU'.i.'.dat' title 'FU'.i,'Files/U'.i.'.dat' title 'U'.i,'Files/FP'.i.'.dat' title 'FP'.i,'Files/P'.i.'.dat' title 'P'.i
# plot 'Files/FU'.i.'.dat' title 'FU'.i,'Files/FP'.i.'.dat' title 'FP'.i
#plot 'Files/FP'.i.'.dat' title 'FP'.i ,'Files/P'.i.'.dat' title 'P'.i

plot 'Files/VP'.i.'.dat' w l title 'VP'.i

i=i+a
if (i < n) reread
