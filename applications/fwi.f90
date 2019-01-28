program main
  use m_kind
  use m_file_function
  use m_polynom
  use m_matrix
  use m_acoustic
  use m_adjoint
  use m_fwi
  use m_init_application
  use m_animation

  implicit none

  !*************** Main Variables ***********************************************
  type(acoustic_problem)              :: forward,backward
  type(t_fwi)                         :: fwi
  real(mp),dimension(:,:),allocatable :: data_P
  real(mp),dimension(:,:),allocatable :: data_U

  real(mp)                            :: t
  integer                             :: i,j
  character(len=50)                   :: file_name

  integer                             :: data_time_step
  integer                             :: fwi_time_step

  integer                             :: values(1:8), k
  integer,dimension(:), allocatable   :: seed
  !******************************************************************************

  call setup_parameters('fwi')

  call date_and_time(values=values)
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)

  !*********** Polynomial initialization *************
  call init_polynom(order)
  
  call init_acoustic_problem(forward,nb_elem,DoF,time_scheme,velocity_data,     &
                             density_data,total_length,final_time,alpha,        &
                             bernstein,signal,boundaries,source_loc,receiver_loc)

  data_time_step=forward%n_time_step
  allocate(data_P(0:forward%n_time_step,2))
  allocate(data_U(0:forward%n_time_step,2))

  t=0
  data_P(0,1)=0.0
  data_P(0,2)=0.0
  if (animation.eq.'data_forward') then
     call print_coef(forward%P,nb_elem,DoF,forward%dx,bernstein,0,'FP')
  end if
  do i=1,forward%n_time_step
     t=i*forward%dt
     call one_forward_step(forward,t)
     data_P(i,1)=t
     data_P(i,2)=forward%P((receiver_loc-1)*DoF+1)
     if (animation.eq.'data_forward')  then
        if (modulo(i,10).eq.0) then
           call print_coef(forward%P,nb_elem,DoF,forward%dx,bernstein,i,'FP')
        end if
     end if
  end do
  data_U=0.0

  open(unit=22,file='data_P.dat',action='write')
  do i=0,forward%n_time_step
     write(22,*) data_P(i,1),data_P(i,2)
  end do
  close(22)

  call free_acoustic_problem(forward)

  print*,' '
  print*,'FWI in progress...'
  print*,' '

  call init_fwi(fwi,nb_iter_fwi,velocity_ini,density_ini,data_P,data_U,      &
                nb_elem,DoF,time_scheme,total_length,final_time,alpha,       &
                bernstein,source_loc,receiver_loc,strategy,scalar_product,   &
                animation,adjoint_test)

  call print_vector(fwi%velocity_model,'VP',0)
  
  do i=1,fwi%nb_iter
     fwi%current_iter=i
     call progress_bar(i,fwi%nb_iter)
     call one_fwi_step(fwi)
     if (animation.eq.'model_update')  then
        if (modulo(i,frame_step).eq.0) then
           call print_vector(fwi%velocity_model,'VP',i)
        end if
     end if
  end do
  call free_fwi(fwi)

  print*,' '
  print*,'... FWI finished'
  print*,' '

  !--------------------- Animation ----------------------------------------------
!  call gif_creation(animation,nb_iter_fwi,data_time_step)
  call gif_creation(gnuplot,animation,frame_step,data_time_step,fwi%nb_iter)
  !------------------------ Free polynoms --------------------------------------
  call free_polynom


end program main
