program main
  use m_file_function
  use m_polynom
  use m_matrix
  use m_acoustic
  use m_adjoint_test

  implicit none

  !*************** Problem Parameters *******************************************
  integer         ,parameter :: nb_elem=200          ! Nb of elements (all same length)
  integer         ,parameter :: ordre=2,DoF=ordre+1  ! Polynoms order
  real            ,parameter :: total_length=1.0     ! domain length
  real            ,parameter :: final_time=0.1       ! final time
  character(len=*),parameter :: time_scheme='LF'    ! change the time scheme
  real            ,parameter :: alpha=1.0            ! Penalisation value
  character(len=*),parameter :: signal='flat'        ! initial values (flat = 0)
  character(len=*),parameter :: boundaries='ABC'     ! Boundary Conditions
  logical         ,parameter :: bernstein=.true.     ! If F-> Lagrange Elements
  integer         ,parameter :: k_max=1e3            ! iter max for power method algo.
  real            ,parameter :: epsilon=1e-5         ! precision for power method algo.
  integer         ,parameter :: n_frame=200         ! nb of times where sol. is saved
  integer         ,parameter :: source_loc=1         ! location of the source (elemts)
  integer         ,parameter :: receiver_loc=2       ! location of the receiver(elemts)
  
  real,dimension(1)          :: velocity ! velocity model change the size to change the model 
  real,dimension(1)          :: density  ! density model change the size to change the model
  !******************************************************************************
  

  !**************** Animation and Outputs ***************************************
  logical,parameter          :: animation=.false.
  logical,parameter          :: sortie=.true. ! animation an RTM not working if F
  logical,parameter          :: RTM=.false.    ! if F -> just forward
  logical,parameter          :: use_data_model=.true.! if T, data = forward receiver
  !******************************************************************************


  !*************** Main Variables ***********************************************
  type(acoustic_problem)        :: forward,backward
  real,dimension(:),allocatable :: P,B,Im,Im_lap   ! vector for Imaging Condition

  integer                           :: i,j
  integer                           :: values(1:8), k
  integer,dimension(:), allocatable :: seed
  real                              :: t0,t1,t2
  logical                           :: file_exists
  !******************************************************************************


  
  !************* Adjoint test****************************************************
  type(t_adjoint_test) :: test
  type(sparse_matrix)  :: Ap
  type(sparse_matrix)  :: Av
  type(sparse_matrix)  :: App
  real                 :: dt
  real                 :: dx
  integer              :: n_adjoint_time_step

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

  
  call cpu_time(t0)
  
  !------------------------------ Forward ---------------------------------------
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%% FORWARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  call init_problem(forward,nb_elem,DoF,time_scheme,velocity,density,            &
                    total_length,final_time,alpha,bernstein,signal,boundaries,  &
                    k_max,epsilon,source_loc,receiver_loc,n_frame,.true.)

  Ap=forward%Ap
  Av=forward%Av
  App=forward%App
  dt=forward%dt
  dx=forward%dx
  n_adjoint_time_step=forward%n_time_step
  
  call print_sol(forward,0)
  call all_time_step(forward,sortie)

  
  if (use_data_model) then
     INQUIRE(FILE="data.dat", EXIST=file_exists)
     if (file_exists) then
        call system('rm data.dat')
        call system('cp receiver.dat data.dat')
     else
        call system('cp receiver.dat data.dat')
     end if
  end if
  
  call cpu_time(t1)
  
  if (RTM) then
     !------------------------------ Backward ------------------------------------
     print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     print*,'%%%%%%%%%%%%%%%%%%%%% BACKWARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

     call init_problem(backward,nb_elem,DoF,time_scheme,velocity,density,        &
          total_length,final_time,alpha,bernstein,signal,boundaries,            &
          k_max,epsilon,source_loc,receiver_loc,n_frame,.false.)
     call print_sol(backward,backward%n_frame)
     call all_time_step(backward,sortie)
  end if
  call cpu_time(t2)


  !---------------------- RTM Post Process ---------------------------------------
  if (RTM) then
     allocate(P(nb_elem*(DoF-1)+1))         ! forward pressure vector
     allocate(B(nb_elem*(DoF-1)+1))         ! backward pressure vector
     allocate(Im(nb_elem*(DoF-1)+1))        ! imaging condition
     allocate(Im_lap(nb_elem*(DoF-1)+1))    ! filtered imaging condition

     Im=0.0

     do i=0,n_frame
        call load(P,'P',i)
        call load(B,'B',i)
        do k=1,nb_elem*(DoF-1)+1
           Im(k)=Im(k)+P(k)*B(k)
        end do
        !     call print_vect(Im,nb_elem,DoF,total_length/nb_elem,.false.,i,'I')
     end do

     Im_lap=laplace_filter(Im)
     open(unit=33,file='RTM.dat')
     open(unit=34,file='RTM_Lap.dat')
     do i=1,nb_elem*(DoF-1)+1
        write(33,*) real(i)/(nb_elem*(DoF-1)+1),Im(i)
        write(34,*) real(i)/(nb_elem*(DoF-1)+1),Im_lap(i)
     end do

     deallocate(P)
     deallocate(B)
     deallocate(Im)
     deallocate(Im_lap)

     call system('gnuplot RTM.script')
     call system('eog RTM.png &')
     !     call system('eog RTM_Lap.png &')
  end if

  
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%% ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

  
  call init_adjoint_test(test,n_adjoint_time_step,time_scheme,dt,Ap,Av,App,nb_elem,DoF,dx)
  call forward_test(test)
  call backward_test(test)
  print*,'inner product UP/DUDP',inner_product(test%P,test%U,test%DP,test%DU)
  print*,'inner product QPQU/FUFP',inner_product(test%QP,test%QU,test%FP,test%FU)
  



  
  print*,'%%%%%%%%%%% END PROGRAM %%%%%%%%%%%%%%%%%%%%%%'
  print*,'Total time = ',t2-t0,'s'
  print*,'Forward time = ',t1-t0,'s'
  print*,'Backward time = ',t2-t1,'s'



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
  call free_acoustic_problem(forward)
  if (RTM) then
     call free_acoustic_problem(backward)
  end if
  call free_basis_b
  call free_basis_l
  call free_B2L
  call free_L2B
  call free_derive



  
end program main
