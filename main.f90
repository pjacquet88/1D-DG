program main
  use load_function
  use m_polynom
  use m_acoustic

  implicit none

  integer:: i,j
  real              :: x,ddx

  integer         ,parameter :: nb_elem=200
  integer         ,parameter :: ordre=2,DoF=ordre+1  ! Polynoms order
  integer         ,parameter :: time_order=2
  real            ,parameter :: total_length=2.0
  real            ,parameter :: final_time=5.0
  real            ,parameter :: alpha=1.0            ! Penalisation value
  character(len=*),parameter :: signal='flat'
  character(len=*),parameter :: boundaries='ABC'
  logical         ,parameter :: bernstein=.true.     ! If F-> Lagrange Elements
  integer         ,parameter :: k_max=1e3
  real            ,parameter :: epsilon=1e-5
  integer         ,parameter :: n_frame=1000
  integer         ,parameter :: receiver_loc=2
  logical         ,parameter :: use_data_model=.true.
  
  real,dimension(1)          :: velocity
  real,dimension(4)          :: density

  
  type(acoustic_problem)     :: problem,forward,backward
  real                       :: t
  integer                    :: n_time_step
  integer                    :: n_display
  real                       :: errorU
  real                       :: errorP
  real,dimension(:),allocatable :: P,B,Im,PB
  integer                       :: nb_frame

  !********************** Animation et Sorties ******************************
  logical,parameter          :: animation=.false.
  logical,parameter          :: sortie=.true.
  !**************************************************************************

  integer :: values(1:8), k
  integer,dimension(:), allocatable :: seed

  real    :: t0,t1,t2

  call date_and_time(values=values)

  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)

  !******************************************************************************
  !********************* ACOUSTIC EQUATION SIMULATION  **************************


  !------------------------ Initialization --------------------------------------
 
  call init_basis_b(ordre)
  call init_basis_l(ordre)

  call create_B2L
  call create_L2B

  call create_derive(ordre)

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
  call init_problem(forward,nb_elem,DoF,time_order,velocity,density,            &
       total_length,final_time,alpha,bernstein,signal,boundaries,  &
                    k_max,epsilon,receiver_loc,n_frame,.true.)
  call print_sol(forward,0)
  call all_time_step(forward,sortie)

  if (use_data_model) then
     call system('rm data.dat')
     call system('cp receiver.dat data.dat')
  end if



  
  call cpu_time(t1)
    !------------------------------ Backward ------------------------------------
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%% BACKWARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  
  call init_problem(backward,nb_elem,DoF,time_order,velocity,density,           &
       total_length,final_time,alpha,bernstein,signal,boundaries,  &
       k_max,epsilon,receiver_loc,n_frame,.false.)

  call print_sol(backward,backward%n_frame)
  
  call all_time_step(backward,sortie)
  call cpu_time(t2)

  
  !--------------------- Animation ----------------------------------------------
  nb_frame=int(forward%n_time_step/forward%n_display)
    if (animation.and.sortie) then
     open(unit=78,file='script.gnuplot',action='write')
     write(78,*)'load "trace1.gnuplot"'
     write(78,*)'n=',nb_frame
     write(78,*)'a=',5
     write(78,*)'load "trace2.gnuplot"'
     close(78)
  
     call system('gnuplot script.gnuplot')
     call system ('eog animate.gif &')
  end if
  

  call error_periodique(problem,t,errorU,errorP)

  !------------------------ Free Variables --------------------------------------
  call free_problem(forward)
  call free_problem(backward)
  call free_basis_b
  call free_basis_l
  call free_B2L
  call free_L2B
  call free_derive

  !---------------------- RTM Post Process ---------------------------------------

  allocate(P(nb_elem*DoF))
  allocate(B(nb_elem*DoF))
  allocate(Im(nb_elem*DoF))

  Im=0.0

  do i=0,n_frame
     call load(P,'P',i)
     call load(B,'B',i)
     do k=1,nb_elem*DoF
        Im(k)=Im(k)+P(k)*B(k)
     end do
!     call print_vect(Im,nb_elem,DoF,total_length/nb_elem,.false.,i,'I')
  end do

  open(unit=33,file='RTM.dat')
  do i=1,nb_elem*DoF
     write(33,*) real(i)/(nb_elem*DoF),Im(i)
  end do

  deallocate(P)
  deallocate(B)
  deallocate(Im)
  
  call system('gnuplot RTM.script')
  call system('eog RTM.png &')
  print*,'%%%%%%%%%%% END RTM %%%%%%%%%%%%%%%%%%%%%%'
  print*,'Total time = ',t2-t0,'s'
  print*,'Forward time = ',t1-t0,'s'
  print*,'Backward time = ',t2-t1,'s'
end program main
