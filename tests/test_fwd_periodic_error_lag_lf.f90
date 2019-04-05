program main
  use m_kind
  use m_file_function
  use m_polynom
  use m_matrix
  use m_acoustic
  use m_init_application
  use m_animation

  implicit none
  
  !*************** Main Variables ***********************************************
  type(acoustic_problem)              :: forward

  real(mp)                            :: t
  integer                             :: i,j,k
  !******************************************************************************

  call setup_parameters('forward')
  bernstein=.false.
  time_scheme='LF'
  signal='ricker'
  boundaries='periodic'
  final_time=1.0_mp
  call print_selected_parameters('forward')
  
  do order=1,10
     DoF=order+1
     call init_polynom(order)
     call init_acoustic_problem(forward,nb_elem,DoF,time_scheme,velocity_data,     &
                                density_data,total_length,final_time,alpha,        &
                                bernstein,flux,signal,boundaries,source_loc,       &
                                receiver_loc)

     t=0

     do i=1,forward%n_time_step
        t=i*forward%dt
        call progress_bar(i,forward%n_time_step)
        call one_forward_step(forward,t)
     end do

     
     print*,'final_time test :',final_time,order,t
  
     
     call free_acoustic_problem(forward)
     call free_polynom
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


! # TEST9 : Forward Periodic error LF
! set(src_fwd_periodic_error_lag_lf tests/test_fwd_periodic_error_lag_lf.f90)
! add_executable(test_fwd_periodic_error_lag_lf ${src_fwd_periodic_error_lag_lf})
! target_link_libraries(test_fwd_periodic_error_lag_lf fwi_lib)
! if(LAPACK_FOUND AND BLAS_FOUND)
!   target_link_libraries(test_fwd_periodic_error_lag_lf ${BLAS_LIBRARIES})
!   target_link_libraries(test_fwd_periodic_error_lag_lf ${LAPACK_LIBRARIES})
! endif()
! add_test(test_fwd_periodic_error_lag_lf ${CMAKE_BINARY_DIR}/test_fwd_periodic_error_lag_lf)
