program main!
  use m_file_function
  use m_polynom
  use m_matrix
  use m_acoustic
  use m_adjoint
  use m_adjoint_test
  use m_fwi

  implicit none

  !*************** Problem Parameters *******************************************
  integer         ,parameter :: nb_elem=100         ! Nb of elements (all same length)
  integer         ,parameter :: ordre=2,DoF=ordre+1  ! Polynoms order
  real            ,parameter :: total_length=1.0     ! domain length
  real            ,parameter :: final_time=3.0      ! final time
  character(len=*),parameter :: time_scheme='RK4'    ! change the time scheme
  real            ,parameter :: alpha=1.0            ! Penalisation value
  character(len=*),parameter :: signal='flat'        ! initial values (flat = 0)
  character(len=*),parameter :: boundaries='ABC'     ! Boundary Conditions
  logical         ,parameter :: bernstein=.false.    ! If F-> Lagrange Elements
  integer         ,parameter :: k_max=1e3            ! itlsr max for power method algo.
  real            ,parameter :: epsilon=1e-5         ! precision for power method algo.
  integer         ,parameter :: n_frame=200          ! nb of times where sol. is saved
  integer         ,parameter :: source_loc=1         ! location of the source (elemts)
  integer         ,parameter :: receiver_loc=30    ! location of the receiver(elemts)

  character(len=*),parameter :: strategy='DTA'
  logical         ,parameter :: adjoint_test=.true.
  character(len=*),parameter :: scalar_product='M'
  integer         ,parameter :: nb_iter_fwi=0
  
  real,dimension(3)          :: velocity ! velocity model change the size to change the model 
  real,dimension(1)          :: density  ! density model change the size to change the model
  !******************************************************************************
  

  !**************** Animation and Outputs ***************************************
  logical,parameter          :: animation=.true.
  logical,parameter          :: sortie=.false. ! animation an RTM not working if F
  logical,parameter          :: RTM=.false.    ! if F -> just forward
  logical,parameter          :: use_data_model=.true.! if T, data = forward receiver
  !******************************************************************************


  !*************** Main Variables ***********************************************
  type(acoustic_problem)            :: forward,backward
  type(t_fwi)                       :: fwi
  real,dimension(:),allocatable     :: P,B,Im,Im_lap   ! vector for Imaging Condition
  real,dimension(:,:),allocatable   :: data_P
  real,dimension(:,:),allocatable   :: data_U

  real                              :: t
  integer                           :: i,j
  integer                           :: values(1:8), k
  integer,dimension(:), allocatable :: seed
  real                              :: t0,t1,t2,t3
  logical                           :: file_exists
  !******************************************************************************

  call date_and_time(values=values)
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)
  !******************************************************************************

  
  !******************************************************************************
  !********************* ACOUSTIC EQUATION SIMULATION  **************************


  !********************** Initialization ***************************************
  
  !*********** Polynomial initialization *************
  call init_basis_b(ordre)
  call init_basis_l(ordre)
  call create_B2L
  call create_L2B
  call create_derive(ordre)
  
  ! do j=1,DoF
  !    do i=0,nb_elem
  !       write(100+j,*) real(i)/real(nb_elem),eval_polynom_l(base_l(j,:),real(i)/real(nb_elem))
  !    end do
  ! end do
  ! STOP
  !***************************************************

  
  !********** Model Initialization *******************
  ! do i=1,size(velocity)
  !    velocity(i)=i
  ! end do
  velocity(1)=1.0
  velocity(2)=1.2
  velocity(3)=1.2
  do i=1,size(density)
     density(i)=1
  end do


  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%% DATA CREATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  
  call init_acoustic_problem(forward,nb_elem,DoF,time_scheme,velocity,density,           &
                    total_length,final_time,alpha,bernstein,signal,boundaries,           &
                    k_max,epsilon,source_loc,receiver_loc)

  allocate(data_P(0:forward%n_time_step,2))
  allocate(data_U(0:forward%n_time_step,2))
  t=0
  data_P(0,1)=0.0
  data_P(0,2)=0.0
!  call print_vect(forward%P,nb_elem,DoF,forward%dx,bernstein,0,'FP')
  do i=1,forward%n_time_step
     t=i*forward%dt
     call one_forward_step(forward,t)
     data_P(i,1)=t
     data_P(i,2)=forward%P((receiver_loc-1)*DoF+2)
 !    call print_vect(forward%P,nb_elem,DoF,forward%dx,bernstein,i,'FP')
  end do
  data_U=0.0

  open(unit=22,file='data_P.dat',action='write')
  do i=0,forward%n_time_step
     write(22,*) data_P(i,1),data_P(i,2)
  end do
  close(22)
  
  call free_acoustic_problem(forward)
 
  call cpu_time(t0)
  

  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%% FWI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

  call init_fwi(fwi,nb_iter_fwi,velocity,density,data_P,data_U,nb_elem,DoF,time_scheme,   &
       total_length,final_time,alpha,bernstein,k_max,epsilon,source_loc,        &
       receiver_loc,strategy,scalar_product,animation,adjoint_test)

  do i=0,fwi%nb_iter

     
    !      do j=1,33
    !    fwi%velocity_model(j)=1.0
    !    fwi%density_model(j)=1.0
    ! end do


    ! do j=34,66
    !    fwi%velocity_model(j)=1.2
    !    fwi%density_model(j)=1.2
    ! end do
    
    ! do j=66,nb_elem
    !    fwi%velocity_model(j)=1.4
    !    fwi%density_model(j)=1.4
    ! end do
    !  if (i.ne.0) then
    !     fwi%velocity_model(i)=fwi%velocity_model(i)+1.0e-5
    !  end if


     
     fwi%current_iter=i
     print*,'Iteration n°',fwi%current_iter,'/',fwi%nb_iter
     call one_fwi_step(fwi)
  end do

  call free_fwi(fwi)



    !--------------------- Animation ----------------------------------------------
  !n_frame=int(forward%n_time_step/forward%n_display)

     ! open(unit=78,file='script.gnuplot',action='write')
     ! write(78,*)'load "trace1.gnuplot"'
     ! write(78,*)'n=',nb_iter_fwi
     ! write(78,*)'a=',1
     ! write(78,*)'load "trace2.gnuplot"'
     ! close(78)
  
     ! call system('gnuplot script.gnuplot')
!     call system ('nw animate.gif &')
    
  !------------------------ Free Variables --------------------------------------
  call free_basis_b
  call free_basis_l
  call free_B2L
  call free_L2B
  call free_derive
  
end program main
