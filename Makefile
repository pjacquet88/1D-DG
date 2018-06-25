FC = gfortran

#CFLAGS = -g -O0 -fbounds-check -ffree-line-length-none -cpp -DKIND_MAT=8 -fdefault-real-8 -freal-4-real-8

CFLAGS = -g -O0 -fbounds-check -ffree-line-length-none -cpp -DKIND_MAT=8

LDFLAG = -llapack -lblas

all: run

%.o: %.f90
	$(FC) $(CFLAGS) -c $<
#	mv $@ .obj

run:   m_matrix.o load_function.o  m_powermethod.o m_polynom.o m_acoustic.o main.o
	$(FC) $(CFLAGS) $^ -o $@  $(LDFLAG)

clean:
	rm -f *.o *.mod *~ *.png fort.* fichier/* animate.gif receiver.dat run
