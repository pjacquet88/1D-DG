program main
  use m_polynom
  use m_acoustic

  implicit none

  integer           :: i,j
  real              :: x,ddx
  type(t_polynom_b) :: bpol
  type(t_polynom_b) :: dbpol
  type(t_polynom_l) :: lpol

  real,parameter      :: h=1.0/10000.0

  integer         ,parameter :: nb_elem=100
  integer         ,parameter :: ordre=4,DoF=ordre+1
  real            ,parameter :: total_length=1.0
  real            ,parameter :: final_time=1.0
  character(len=*),parameter :: signal='sinus'
  character(len=*),parameter :: boundaries='periodique'
  logical         ,parameter :: bernstein=.true.
  type(acoustic_problem)     :: problem

  integer                    :: n_time_step
  integer         ,parameter :: n_display=10

  real,dimension(DoF) :: coef_test_b,coef_test_l
  real,dimension(DoF*nb_elem) :: work

  integer :: values(1:8), k
  integer, dimension(:), allocatable :: seed

  call date_and_time(values=values)

  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)
  

  call init_basis_b(ordre)
  call init_basis_l(ordre)

  call create_B2L
  call create_L2B



  !****************SORTIE VISUELLE DES FONCTIONS DE BASES***********************
  do i=0,100
     do j=1,ordre+1
        write(100+j,*),i/100.0,eval_polynom_b(base_b(j),i/100.0)
        write(200+j,*),i/100.0,eval_polynom_l(base_l(j),i/100.0)
     end do
  end do
  call system('gnuplot script_base_b.gnu')
  call system('gnuplot script_base_l.gnu')
  !*****************************************************************************


  !******************** TEST B2L ***********************************************

  do i=1,ordre+1
     call random_number(coef_test_b(i))
  end do

  
  call init_polynom_b(bpol,coef_test_b)

  call Bernstein2Lagrange(bpol,lpol)
  
  do i=0,100
     write(20,*),i/100.0,eval_polynom_b(bpol,i/100.0)
     write(21,*),i/100.0,eval_polynom_l(lpol,i/100.0)
  end do

  call system('gnuplot script_test_b2l.gnu')
  !call system('eog b2l.png &')

  call free_polynom_b(bpol)
  call free_polynom_l(lpol)

  !******************** TEST L2B ***********************************************

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
 ! call system('eog l2b.png &')

  call free_polynom_b(bpol)
  call free_polynom_l(lpol)

!****************** TEST DERIVE *******************************************

  do i=1,ordre+1
     call random_number(coef_test_b(i))
  end do


  call init_polynom_b(bpol,coef_test_b)

  call create_derive(ordre)
  call deriv_pol_b(bpol,dbpol)


   write(31,*) 0.0, (eval_polynom_b(bpol,0.0+h)-eval_polynom_b(bpol,0.0))/(h)
     write(30,*) 0.0, eval_polynom_b(dbpol,0.0)

  do i=1,99
     write(31,*) i/100.0, (eval_polynom_b(bpol,i/100.0+h)-eval_polynom_b(bpol,i/100.0-h))/(2*h)
     write(30,*) i/100.0, eval_polynom_b(dbpol,i/100.0)
  end do

     write(31,*) 1.0, (eval_polynom_b(bpol,1.0)-eval_polynom_b(bpol,1.0-h))/(h)
     write(30,*) 1.0, eval_polynom_b(dbpol,1.0)


  call system('gnuplot script_test_derive.gnu')
  !call system('eog derive.png &')

  !*************************************************************************
  !********************* SIMULATION EQUATION ACOUSTIQUE*********************

  
  call init_quadrature
  call init_problem(problem,nb_elem,DoF,total_length,final_time,bernstein,signal)
  call print_sol(problem,0)


  call init_m_loc(m_loc,m_inv_loc,DoF,bernstein)
  call init_s_loc(s_loc,DoF,bernstein)
  call init_stiffness(Av,Ap,m_inv_loc,s_loc,nb_elem,DoF,boundaries,bernstein,problem)

  print*,'le pas d espace est de :',problem%dx
  n_time_step=int(final_time/problem%dt)
  print*,'il y aura :',n_time_step,'time _step'

  call one_time_step(problem)
  call print_sol(problem,1)

  ! do i=1,n_time_step
  !    call one_time_step(problem)

  !    if (modulo(i,n_display).eq.0) then
  !       call print_sol(problem,i/n_display)
  !    end if
  ! end do


  ! erreur
  ! ddx=problem%dx/(DoF-1)
  ! open (unit=4,file='fichier/error.dat',action='write')
  ! do i=1,nb_elem
  !    x=(i-1)*problem%dx-ddx
  !    do j=1,DoF
  !       x=x+ddx
  !       work((i-1)*DoF+j)=signal_ini(x-final_time,0.0,signal)
  !       write(4,*) x,work((i-1)*DoF+j)
  !    end do
  ! end do

  ! close(4)

  ! call system('gnuplot test.gnu')
  ! call system('eog test.png &')

  ! open(unit=78,file='script.gnuplot',action='write')
  ! write(78,*),'load "trace1.gnuplot"'
  ! write(78,*),'n=',n_time_step
  ! write(78,*),'a=',1
  ! write(78,*),'load "trace2.gnuplot"'
  ! close(78)

  ! call system('gnuplot5 script.gnuplot')
  
  print*,'done'
  
end program main
