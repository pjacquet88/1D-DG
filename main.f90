program main
  use m_polynom

  implicit none


  integer           :: i,j
  type(t_polynom_b) :: bpol
  type(t_polynom_b) :: bpol1,bpol2,bpol3,bpol4,bpol5
  type(t_polynom_b) :: dbpol
  type(t_polynom_l) :: lpol
  type(t_polynom_l) :: lpol1,lpol2,lpol3,lpol4,lpol5

  real,dimension(5) :: coef_test_b
  real,dimension(5) :: coef_test_db
  real,dimension(5) :: coef_test_l
  real              :: erreur,a

  real,dimension(5),parameter :: coef1=(/1.0,0.0,0.0,0.0,0.0/)
  real,dimension(5),parameter :: coef2=(/0.0,1.0,0.0,0.0,0.0/)
  real,dimension(5),parameter :: coef3=(/0.0,0.0,1.0,0.0,0.0/)
  real,dimension(5),parameter :: coef4=(/0.0,0.0,0.0,1.0,0.0/)
  real,dimension(5),parameter :: coef5=(/0.0,0.0,0.0,0.0,1.0/)


  real,parameter :: h=1/10000.0
    
  integer :: values(1:8), k
  integer, dimension(:), allocatable :: seed

  call date_and_time(values=values)

  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)

  print*,'test0'

  call init_polynom_b(bpol1,coef1)
  call init_polynom_b(bpol2,coef2)
  call init_polynom_b(bpol3,coef3)
  call init_polynom_b(bpol4,coef4)
  call init_polynom_b(bpol5,coef5)

  call init_polynom_l(lpol1,coef1)
  call init_polynom_l(lpol2,coef2)
  call init_polynom_l(lpol3,coef3)
  call init_polynom_l(lpol4,coef4)
  call init_polynom_l(lpol5,coef5)

  print*,'test1'
  
  do i=0,100
     write(21,*),i/100.0,eval_polynom_b(bpol1,i/100.0)
     write(22,*),i/100.0,eval_polynom_b(bpol2,i/100.0)
     write(23,*),i/100.0,eval_polynom_b(bpol3,i/100.0)
     write(24,*),i/100.0,eval_polynom_b(bpol4,i/100.0)
     write(25,*),i/100.0,eval_polynom_b(bpol5,i/100.0)

     write(31,*),i/100.0,eval_polynom_l(lpol1,i/100.0)
     write(32,*),i/100.0,eval_polynom_l(lpol2,i/100.0)
     write(33,*),i/100.0,eval_polynom_l(lpol3,i/100.0)
     write(34,*),i/100.0,eval_polynom_l(lpol4,i/100.0)
     write(35,*),i/100.0,eval_polynom_l(lpol5,i/100.0)
  end do

  call system('gnuplot script_base_b.gnu')
  !call system('eog base_b.png &')

  call system('gnuplot script_base_l.gnu')
  !call system('eog base_l.png &')
  
  
  do i=1,5
     call random_number(coef_test_b(i))
  end do

  call init_polynom_b(bpol,coef_test_b)

  print*,'test2'
  call create_B2L
  call Bernstein2Lagrange(bpol,lpol)

  do i=0,100
     write(41,*),i/100.0,eval_polynom_b(bpol,i/100.0)
     write(42,*),i/100.0,eval_polynom_l(lpol,i/100.0)
  end do

  call system('gnuplot script_test_b2l.gnu')
 ! call system('eog b2l.png &')

  do i=0,100

  end do

  call create_derive
  call deriv_pol_b(bpol,dbpol)

  print*,'test3'
  
  do i=0,100
     write(51,*) i/100.0, (eval_polynom_b(bpol,i/100.0+h)-eval_polynom_b(bpol,i/100.0-h))/(2*h)
     write(52,*) i/100.0, eval_polynom_b(dbpol,i/100.0)
  end do


 call system('gnuplot script_test_derive.gnu')
 call system('eog derive.png &')

 call create_L2B
 
  print*,'done'

  
end program main
