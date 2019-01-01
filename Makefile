FC = gfortran


# All real in real*8
CFLAGS = -g -O3 -fbounds-check -ffree-line-length-none -cpp -DKIND_MAT=8 -fdefault-real-8 -freal-4-real-8


# All real in real*4
#CFLAGS = -g -O0 -fbounds-check -ffree-line-length-none -cpp -DKIND_MAT=8

LDFLAG = -llapack -lblas

all: run

%.o: %.f90
	$(FC) $(CFLAGS) -c $<
#	mv $@ .obj

run:    m_matrix.o m_time_scheme.o m_powermethod.o m_polynom.o m_file_function.o m_adjoint_test.o m_acoustic.o m_adjoint.o m_fwi.o  m_init_application.o m_animation.o main.o
	$(FC) $(CFLAGS) $^ -o $@  $(LDFLAG)

clean:
	rm -f *.o *.mod *~ *.png fort.* Files/* animate.gif *.dat run
