#plot 'fichier/P'.i.'.dat' title 'P'.i
#plot 'fichier/B'.i.'.dat' title 'B'.i
plot 'fichier/P'.i.'.dat' title 'P'.i,'fichier/B'.i.'.dat' title 'B'.i
i=i+a
if (i < n) reread
