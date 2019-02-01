program main
  use m_kind
  use m_polynom
  use m_matrix
  use m_init_application
  use m_acoustic

  implicit none

  !*************** Main Variables ***********************************************
  type(acoustic_problem)              :: forward
  integer                             :: i,j,k
  real(mp),dimension(:),allocatable   :: x_lag,x_bb
  logical                             :: bool
  real(mp)                            :: error
  !******************************************************************************
 
  call setup_parameters('forward')
  bool=.true.
  call print_selected_parameters('forward')
  
  do order=1,10
     DoF=order+1
     call init_polynom(order)
     allocate(x_lag(DoF*nb_elem))
     allocate(x_bb(DoF*nb_elem))
     x_lag=1.0_mp
     
     do i=1,nb_elem
        call Lagrange2Bernstein(x_lag((i-1)*DoF+1:i*DoF),x_bb((i-1)*DoF+1:i*DoF))
     end do

     call init_acoustic_problem(forward,nb_elem,DoF,time_scheme,velocity_data,     &
          density_data,total_length,final_time,alpha,        &
          bernstein,signal,boundaries,source_loc,receiver_loc)

     error=abs(dot_product(sparse_matmul(forward%S,x_bb),x_bb)/nb_elem)

     print*,'order ::',order,'error ::',error

     if (mp.eq.sp) then
        bool=(bool.and.(error.lt.1e-4_mp))
     else if (mp.eq.dp) then
        bool=(bool.and.(error.lt.1e-8_mp))
     end if

     call my_return(bool)
     call free_acoustic_problem(forward)
     call free_polynom
     deallocate(x_lag)
     deallocate(x_bb)
  end do

contains

  subroutine my_return(bool)
    logical,intent(in) :: bool

    if (bool) then
       print*,'Succeed'
    else
       ERROR STOP "TEST FAILED"
    end if
  end subroutine my_return

end program main
