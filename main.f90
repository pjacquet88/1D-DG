program main
  use m_polynom
  use m_acoustic

  implicit none

  integer           :: i,j
  real              :: x,ddx

  !  integer         ,parameter :: nb_elem=100
  integer                    :: nb_elem
  integer         ,parameter :: ordre=4,DoF=ordre+1  ! Polynoms order
  !integer         ,parameter :: ordre,DoF
  real            ,parameter :: total_length=1.0
  real            ,parameter :: final_time=1.0
  real            ,parameter :: alpha=0.5            ! Penalisation value
  character(len=*),parameter :: signal='ricker'
  character(len=*),parameter :: boundaries='periodique'
  logical         ,parameter :: bernstein=.true.    !If F-> Lagrange Elements
  logical         ,parameter :: F_forte=.true. !Means Strong Formulation (always=T)
  integer         ,parameter :: k_max=1e3
  real            ,parameter :: epsilon=1e-5
  type(acoustic_problem)     :: problem
  
  real                       :: t
  integer                    :: n_time_step
  integer                    :: n_display
  real                       :: errorU
  real                       :: errorP
  real                       :: test

  !********************** Animation et Sorties ******************************
  logical,parameter          :: animation=.false.
  logical,parameter          :: sortie=.false.
  !**************************************************************************

  integer :: values(1:8), k
  integer,dimension(:), allocatable :: seed

  call date_and_time(values=values)

  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)

  !*************************************************************************
  !********************* SIMULATION EQUATION ACOUSTIQUE*********************



  ! do ordre=1,10
  !    DoF=ordre+1

  call init_basis_b(ordre)
  call init_basis_l(ordre)

  print*,'test'
  do i=1,ordre+1
     print*,base_l(i,:)
  end do

  call create_B2L
  call create_L2B


  call create_derive(ordre)

  do nb_elem=10,100,1

     call init_problem(problem,nb_elem,DoF,total_length,final_time,alpha,bernstein,&
          F_forte,signal,boundaries,k_max,epsilon)

     call init_ApAv(problem)
     call init_UP(problem)
     call print_sol(problem,0)
     !call error_periodique(problem,t,errorU,errorP,0)
     n_time_step=int(final_time/problem%dt)
     if (n_time_step.le.200) then
        n_display=1
     else
        n_display=n_time_step/200
     end if
     print*,'n_display',n_display
     print*,'dt = ',problem%dt
     print*,'There will be :',n_time_step,'time_step'

     call one_time_step(problem,t)
     call print_sol(problem,1)
     
     do i=2,n_time_step
        t=i*problem%dt
        call one_time_step(problem,t)

        if (sortie) then
           if (modulo(i,n_display).eq.0) then
              call print_sol(problem,i/n_display)
              !call error_periodique(problem,t,errorU,errorP,i/n_display)
             ! write(88,*) t,log10(errorU)
           end if
        end if
           ! if (i.eq.n_time_step) then
           !    call print_sol(problem,42)
           ! end if

        end do

        if (animation.and.sortie) then
           open(unit=78,file='script.gnuplot',action='write')
           write(78,*)'load "trace1.gnuplot"'
           write(78,*)'n=',int(n_time_step/n_display)
           write(78,*)'a=',1
           write(78,*)'load "trace2.gnuplot"'
           close(78)

           call system('gnuplot script.gnuplot')
           call system ('eog animate.gif &')
        end if


        call error_periodique(problem,t,errorU,errorP)
        call free_problem(problem)
        write(88,*) log10(1.0/nb_elem),log10(errorU)
        write(89,*) log10(1.0/nb_elem),log10(errorP)
        ! write(77,*) ordre,problem%dt

     end do


     call free_basis_b
     call free_basis_l

     call free_B2L
     call free_L2B
     call free_derive

   end program main
