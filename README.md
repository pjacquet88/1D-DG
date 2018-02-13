# Bernstein

To compile : "make"
To clean : "make clean"
To execute : "./run"

You need to: -create an empty folder called 'fichier' at the source code root
             -have a recent version of gnuplot to create and see the animation automatically created
             -check your Blas/Lapack library location in the Makefile
            
            If there is a bug with gnuplot replace 'gnuplot5' at the end of the main.f90 by 'gnuplot'
            
 Some parameters at the beginning of the file main.f90 can changed in order to personalize the simulation

There is a excel file (comparaison.xls) that compare the "sparcity" of the stiffness matrix of Lagrange and Bernstain polynoms.
 
 NB : no works has been done for the cfl. That is to say that the step time length dt is not actualised with the degree of polynom.
 you need to change directly the value in the file m_acoustic.f90 (search for "dt=")
