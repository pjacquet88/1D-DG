# Bernstein

To compile : "make"
To clean : "make clean"
To execute : "./run"

You need to: -create an empty folder called 'Files' at the source code root ($ mkdir Files)
             -have a recent version of gnuplot to create and see the animation automatically created ($ sudo apt install gnuplot)
             -check your Blas/Lapack library location in the Makefile
            
            
 Some parameters at the beginning of the file main.f90 can be changed in order to personalize the simulation
For instance :
 
_ RTM=T leads to a RTM simulation
_ RTM=F leads to a forward simulation
_ Animation = T or F to have a forward animation (the backward propagation can be obtained by modifying the script animate.gnuplot)
_ the size of the vector density and velocity can be changed manually in order to create the desired amount of interfaces. If the size of velocity is X it will cut the entire domain in X equal part of different velocity (in the main the velocity is given by an incrementation but it can be personnalized however you want)

By default the simulation is perfomed in simple precision. If you want to switch to double precision comment/uncomment the CFLAGS lines and change in the m_matrix.f90 module (line 33/34) SGETRF and SGETRI to DGETRF and DGETRI 
 