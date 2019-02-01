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
  integer                             :: forward_n_time_step

  logical                             :: bool
  !******************************************************************************
  bool=.true.
  
  call setup_parameters('forward')
  final_time=10
  velocity_data=1.0_mp
  time_scheme='AB3'
  bernstein=.true.
  signal='flat'
  animation='no'
  gnuplot=.false.
  call print_selected_parameters('forward')


  ! ----------- BERNSTEIN ------------------------------------
    print*,'TEST ON BERNSTEIN'
  call init_polynom(order)

  call init_acoustic_problem(forward,nb_elem,DoF,time_scheme,velocity_data,     &
                             density_data,total_length,final_time,alpha,        &
                             bernstein,signal,boundaries,source_loc,receiver_loc)

  forward_n_time_step=forward%n_time_step

  t=0
  if (animation.eq.'data_forward') then
     call print_coef(forward%P,nb_elem,DoF,forward%dx,bernstein,0,'FP')
  end if
  do i=1,forward%n_time_step
     t=i*forward%dt
     call progress_bar(i,forward%n_time_step)
     call one_forward_step(forward,t)
     if (animation.eq.'data_forward')  then
        if (modulo(i,frame_step).eq.0) then
           call print_coef(forward%P,nb_elem,DoF,forward%dx,bernstein,i,'FP')
        end if
     end if
  end do

  print*,'test bernstein :'
  print*,'minval/maxval P :',minval(forward%P),maxval(forward%P)
  print*,'minval/maxval U :',minval(forward%U),maxval(forward%U)
  
  bool=bool.and.(minval(abs(forward%P)).le.1e-5)
  bool=bool.and.(maxval(abs(forward%P)).le.1e-5)
  bool=bool.and.(minval(abs(forward%U)).le.1e-5)
  bool=bool.and.(maxval(abs(forward%U)).le.1e-5)

  call my_return(bool)
  
  call free_acoustic_problem(forward)
  call free_polynom

  bernstein=.false.
  ! ----------- LAGRANGE ------------------------------------
  call print_selected_parameters('forward')
  print*,'TEST ON LAGRANGE'
  call init_polynom(order)

  call init_acoustic_problem(forward,nb_elem,DoF,time_scheme,velocity_data,     &
                             density_data,total_length,final_time,alpha,        &
                             bernstein,signal,boundaries,source_loc,receiver_loc)

  forward_n_time_step=forward%n_time_step

  t=0
  if (animation.eq.'data_forward') then
     call print_coef(forward%P,nb_elem,DoF,forward%dx,bernstein,0,'FP')
  end if
  do i=1,forward%n_time_step
     t=i*forward%dt
     call progress_bar(i,forward%n_time_step)
     call one_forward_step(forward,t)
     if (animation.eq.'data_forward')  then
        if (modulo(i,frame_step).eq.0) then
           call print_coef(forward%P,nb_elem,DoF,forward%dx,bernstein,i,'FP')
        end if
     end if
  end do

  print*,'test lagrange :'
  print*,'minval/maxval P :',minval(forward%P),maxval(forward%P)
  print*,'minval/maxval U :',minval(forward%U),maxval(forward%U)
  
  bool=bool.and.(minval(abs(forward%P)).le.1e-5)
  bool=bool.and.(maxval(abs(forward%P)).le.1e-5)
  bool=bool.and.(minval(abs(forward%U)).le.1e-5)
  bool=bool.and.(maxval(abs(forward%U)).le.1e-5)

  call my_return(bool)
  
  call free_acoustic_problem(forward)
  call free_polynom

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
