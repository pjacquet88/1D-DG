program main
  use m_polynom
  use m_acoustic

  implicit none

  integer           :: i,j
  type(t_polynom_b) :: bpol
  type(t_polynom_b) :: dbpol
  type(t_polynom_l) :: lpol

  real,parameter :: h=1/100000.0
  integer,parameter :: ordre=10,DoF=ordre+1

  real,dimension(DoF) :: coef_test_b
  real,dimension(DoF) :: coef_test_db
  real,dimension(DoF) :: coef_test_l
  real                :: erreur,a

  real,dimension(ordre+1,ordre+1) :: matrix
  real,dimension(ordre+1,ordre+1) :: inv_matrix
  real,dimension(ordre+1,ordre+1) :: test_matrix

  integer :: values(1:8), k
  integer, dimension(:), allocatable :: seed

  call date_and_time(values=values)

  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)
  

  call init_basis_b(ordre)
  call init_basis_l(ordre)


  do i=0,100
     do j=1,ordre+1
        write(100+j,*),i/100.0,eval_polynom_b(base_b(j),i/100.0)
        write(200+j,*),i/100.0,eval_polynom_l(base_l(j),i/100.0)
     end do
  end do
  call system('gnuplot script_base_b.gnu')
  call system('gnuplot script_base_l.gnu')


  call create_B2L
  call create_L2B

  do i=1,ordre+1
     call random_number(coef_test_l(i))
  end do



  call init_polynom_l(lpol,coef_test_l)

  call Lagrange2Bernstein(lpol,bpol)

  do i=0,100
     write(10,*),i/100.0,eval_polynom_b(bpol,i/100.0)
     write(11,*),i/100.0,eval_polynom_l(lpol,i/100.0)
  end do

  call system('gnuplot script_test_l2b.gnu')

  call system('eog l2b.png &')

  do i=0,100
     write(22,*),i/100.0,signal_ini(i/100.0,'ricker')
  end do


  print*,'done'

end program main
