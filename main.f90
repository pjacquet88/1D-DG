program main
  use m_polynom
  use m_acoustic

  implicit none

  integer           :: i,j
  type(t_polynom_b) :: bpol
  type(t_polynom_b) :: dbpol
  type(t_polynom_l) :: lpol

  real,parameter :: h=1/100000.0
  integer,parameter :: ordre=12,DoF=ordre+1

  real,dimension(DoF) :: coef_test_b
  real,dimension(DoF) :: coef_test_db
  real,dimension(DoF) :: coef_test_l
  real                :: erreur,a
  
  integer :: values(1:8), k
  integer, dimension(:), allocatable :: seed

  call date_and_time(values=values)

  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)


   call init_basis_b(ordre)
   call init_basis_l(ordre)
   
   
   ! do j=1,ordre+1
   !    do i=0,100
   !       write(20+j,*),i/100.0,eval_polynom_b(base_b(j),i/100.0)
   !       write(30+j,*),i/100.0,eval_polynom_l(base_l(j),i/100.0)
   !    end do
   ! end do

  !call system('gnuplot script_base_b.gnu')
  !call system('eog base_b.png &')

  !call system('gnuplot script_base_l.gnu')
  !call system('eog base_l.png &')


  do i=1,ordre+1
     call random_number(coef_test_b(i))
  end do

  call init_polynom_b(bpol,coef_test_b)

  call create_B2L
  
  call Bernstein2Lagrange(bpol,lpol)

  do i=0,100
     write(141,*),i/100.0,eval_polynom_b(bpol,i/100.0)
     write(142,*),i/100.0,eval_polynom_l(lpol,i/100.0)
  end do

  ! call free_polynom_l(lpol)
  ! call free_polynom_b(bpol)

  call system('gnuplot script_test_b2l.gnu')
  !call system('eog b2l.png &')

  ! do i=1,ordre+1
  !    call random_number(coef_test_l(i))
  ! end do

  ! call init_polynom_l(lpol,coef_test_l)

  ! call create_L2B
  
  ! call Lagrange2Bernstein(lpol,bpol)

  ! do i=0,100
  !    write(151,*),i/100.0,eval_polynom_l(lpol,i/100.0)
  !    write(152,*),i/100.0,eval_polynom_b(bpol,i/100.0)
  ! end do

  ! call system('gnuplot script_test_l2b.gnu')
  ! call system('eog l2b.png &')

  call create_derive(ordre)
  call deriv_pol_b(bpol,dbpol)


   write(251,*) 0.0, (eval_polynom_b(bpol,0.0+h)-eval_polynom_b(bpol,0.0))/(h)
     write(252,*) 0.0, eval_polynom_b(dbpol,0.0)

  do i=1,99
     write(251,*) i/100.0, (eval_polynom_b(bpol,i/100.0+h)-eval_polynom_b(bpol,i/100.0-h))/(2*h)
     write(252,*) i/100.0, eval_polynom_b(dbpol,i/100.0)
  end do

     write(251,*) 1.0, (eval_polynom_b(bpol,1.0)-eval_polynom_b(bpol,1.0-h))/(h)
     write(252,*) 1.0, eval_polynom_b(dbpol,1.0)



 call system('gnuplot script_test_derive.gnu')
 call system('eog derive.png &')

 ! call create_L2B

 ! call init_quadrature
 ! call init_m_loc_l(DoF)

 ! do i=1,DoF
 !    print*,m_loc(i,:)
 ! end do


 print*,'done'


 

end program main
