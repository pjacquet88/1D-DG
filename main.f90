program main!
  use m_file_function
  use m_polynom
  use m_matrix
  use m_acoustic
  use m_adjoint_test
  use m_fwi

  implicit none

  !*************** Problem Parameters *******************************************
  integer         ,parameter :: nb_elem=150          ! Nb of elements (all same length)
  integer         ,parameter :: ordre=2,DoF=ordre+1  ! Polynoms order
  real            ,parameter :: total_length=1.0     ! domain length
  real            ,parameter :: final_time=2.0       ! final time
  character(len=*),parameter :: time_scheme='AB3'    ! change the time scheme
  real            ,parameter :: alpha=1.0            ! Penalisation value
  character(len=*),parameter :: signal='flat'        ! initial values (flat = 0)
  character(len=*),parameter :: boundaries='ABC'     ! Boundary Conditions
  logical         ,parameter :: bernstein=.true.     ! If F-> Lagrange Elements
  integer         ,parameter :: k_max=1e3            ! itlsr max for power method algo.
  real            ,parameter :: epsilon=1e-5         ! precision for power method algo.
  integer         ,parameter :: n_frame=200          ! nb of times where sol. is saved
  integer         ,parameter :: source_loc=1         ! location of the source (elemts)
  integer         ,parameter :: receiver_loc=10      ! location of the receiver(elemts)
  
  real,dimension(1)          :: velocity ! velocity model change the size to change the model 
  real,dimension(3)          :: density  ! density model change the size to change the model
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
  !***************************************************

  
  !********** Model Initialization *******************
  do i=1,size(velocity)
     velocity(i)=i
  end do
  do i=1,size(density)
     density(i)=i
  end do

  call init_problem(forward,nb_elem,DoF,time_scheme,velocity,density,           &
                    total_length,final_time,alpha,bernstein,signal,boundaries,  &
                    k_max,epsilon,source_loc,receiver_loc,1,.true.)

  allocate(data_P(0:forward%n_time_step,2))
  allocate(data_U(0:forward%n_time_step,2))
  t=0
  data_P(0,1)=0.0
  data_P(0,2)=0.0
  do i=1,forward%n_time_step
     t=i*forward%dt
     call one_time_step(forward,t)
     data_P(i,1)=t
     data_P(i,2)=forward%P(receiver_loc)
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
  print*,'%%%%%%%%%%%%%%%%%%%%% FWI %%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

  call init_fwi(fwi,1,velocity,density,data_P,data_U,nb_elem,DoF,time_scheme,   &
       total_length,final_time,alpha,bernstein,k_max,epsilon,source_loc,        &
       receiver_loc)

  call one_fwi_step(fwi)

  call free_fwi(fwi)



    !--------------------- Animation ----------------------------------------------
  !n_frame=int(forward%n_time_step/forward%n_display)
    if (animation.and.sortie) then
     open(unit=78,file='script.gnuplot',action='write')
     write(78,*)'load "trace1.gnuplot"'
     write(78,*)'n=',n_frame
     write(78,*)'a=',5
     write(78,*)'load "trace2.gnuplot"'
     close(78)
  
     call system('gnuplot script.gnuplot')
     call system ('eog animate.gif &')
  end if
    
  !------------------------ Free Variables --------------------------------------
  call free_basis_b
  call free_basis_l
  call free_B2L
  call free_L2B
  call free_derive
  
end program main
