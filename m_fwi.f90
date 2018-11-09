module m_fwi
  use m_file_function
  use m_matrix
  use m_time_scheme
  use m_acoustic
  use m_adjoint
  use m_adjoint_test
  implicit none


  type t_fwi
     type(acoustic_problem)          :: forward
     type(adjoint_problem)           :: backward
     integer                         :: nb_iter
     integer                         :: current_iter
     integer                         :: model_size
     real,dimension(:)  ,allocatable :: velocity_model
     real,dimension(:)  ,allocatable :: density_model
     integer                         :: n_time_step
     real,dimension(:,:),allocatable :: data_P
     real,dimension(:,:),allocatable :: data_U
     real,dimension(:,:),allocatable :: P_received
     real,dimension(:,:),allocatable :: U_received
     real,dimension(:,:),allocatable :: P
     real,dimension(:,:),allocatable :: U
     real,dimension(:,:),allocatable :: QP
     real,dimension(:,:),allocatable :: QU
     real,dimension(:,:),allocatable :: FP
     real,dimension(:,:),allocatable :: FU
     real,dimension(:,:),allocatable :: FP_half
     real,dimension(:,:),allocatable :: FU_half
     real,dimension(:,:),allocatable :: DP
     real,dimension(:,:),allocatable :: DU
     real,dimension(:,:),allocatable :: GP
     real,dimension(:,:),allocatable :: GU
     real,dimension(:)  ,allocatable :: gradient_rho
     real,dimension(:)  ,allocatable :: gradient_vp
     real,dimension(:,:),allocatable :: M
     real,dimension(:,:),allocatable :: Ap
     real,dimension(:,:),allocatable :: App
     real,dimension(:,:),allocatable :: Av
     
     integer                         :: nb_elem
     integer                         :: DoF
     character(len=20)               :: time_scheme
     real                            :: total_length
     real                            :: final_time
     real                            :: dt
     real                            :: dx
     real                            :: alpha
     logical                         :: bernstein
     character(len=20)               :: signal
     character(len=20)               :: boundaries
     integer                         :: k_max
     real                            :: epsilon
     integer                         :: source_loc
     integer                         :: receiver_loc
     character(len=20)               :: strategy
     character(len=20)               :: scalar_product
     logical                         :: adjoint_test

     logical                         :: animation
  end type t_fwi
    

contains

  subroutine init_fwi(fwi,nb_iter,velocity_model,density_model,data_P,data_U,   &
                      nb_elem,DoF,time_scheme,total_length,final_time,alpha,    &
                      bernstein,k_max,epsilon,source_loc,receiver_loc,strategy, &
                      scalar_product,animation,adjoint_test)
    type(t_fwi)        ,intent(inout) :: fwi
    integer            ,intent(in)    :: nb_iter
    real,dimension(:)  ,intent(in)    :: velocity_model
    real,dimension(:)  ,intent(in)    :: density_model
    real,dimension(:,:),intent(in)    :: data_P
    real,dimension(:,:),intent(in)    :: data_U
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    character(len=*)   ,intent(in)    :: time_scheme
    real               ,intent(in)    :: total_length
    real               ,intent(in)    :: final_time
    real               ,intent(in)    :: alpha
    logical            ,intent(in)    :: bernstein
    integer            ,intent(in)    :: k_max
    real               ,intent(in)    :: epsilon
    integer            ,intent(in)    :: source_loc
    integer            ,intent(in)    :: receiver_loc
    character(len=*)   ,intent(in)    :: strategy
    character(len=*)   ,intent(in)    :: scalar_product
    logical            ,intent(in)    :: animation
    logical            ,intent(in)    :: adjoint_test
    
    integer :: i
    real    :: epsilon_vp
    character(len=50) :: fichier
    fwi%nb_iter=nb_iter
    fwi%model_size=size(density_model)

    ! allocate(fwi%velocity_model(fwi%model_size))
    ! allocate(fwi%density_model(fwi%model_size))
    ! fwi%velocity_model=velocity_model
    ! fwi%density_model=density_model

    allocate(fwi%velocity_model(nb_elem))
    allocate(fwi%density_model(nb_elem))

    epsilon_vp=10.0d0**(-5)
    
    do i=1,33
       fwi%velocity_model(i)=1.0!+epsilon_vp
       fwi%density_model(i)=1.0!+10**(-5)
    end do


    do i=34,66
       fwi%velocity_model(i)=1.0!epsilon_vp
       fwi%density_model(i)=1.0!+10**(-5)
    end do
    
    do i=67,nb_elem
       fwi%velocity_model(i)=1.0!+epsilon_vp
       fwi%density_model(i)=1.0!+10**(-5)
    end do

    
    allocate(fwi%gradient_rho(nb_elem))
    allocate(fwi%gradient_vp(nb_elem))

    allocate(fwi%data_P(0:size(data_P,1)-1,size(data_U,2)))
    allocate(fwi%data_U(0:size(data_U,1)-1,size(data_U,2)))

    fwi%data_P=data_P
    fwi%data_U=data_U

    fwi%nb_elem=nb_elem
    fwi%DoF=DoF
    fwi%time_scheme=time_scheme
    fwi%total_length=total_length
    fwi%final_time=final_time
    fwi%alpha=alpha
    fwi%bernstein=bernstein
    fwi%k_max=k_max
    fwi%epsilon=epsilon
    fwi%source_loc=source_loc
    fwi%receiver_loc=receiver_loc
    fwi%strategy=strategy
    fwi%scalar_product=scalar_product
    fwi%animation=animation
    fwi%adjoint_test=adjoint_test
    
    fwi%signal='flat'
    fwi%boundaries='ABC'
    
    write(fichier,"(A,A,I0,'.dat')") "../Bernstein/Files/",'VP',0
    open(unit=28,file=fichier)
    do i=1,fwi%nb_elem
       write(28,*) i,fwi%velocity_model(i)
    end do


  end subroutine init_fwi
  

  subroutine free_fwi(fwi)
    type(t_fwi),intent(inout) :: fwi

    deallocate(fwi%velocity_model)
    deallocate(fwi%density_model)
    deallocate(fwi%data_P)
    deallocate(fwi%data_U)
    deallocate(fwi%gradient_vp)
    deallocate(fwi%gradient_rho)
  end subroutine free_fwi


  subroutine one_fwi_step(fwi)
    type(t_fwi),intent(inout) :: fwi
    integer                   :: current_time_step
    real                      :: t
    real                      :: t1,t2

    real,dimension(fwi%DoF*fwi%nb_elem) :: RTM
    integer                             :: i,j
    type(t_adjoint_test)                :: test

    !------------------------------ Forward -------------------------------------
    ! print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    ! print*,'%%%%%%%%%%%%%%%%%%%%% FORWARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    ! print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call init_acoustic_problem(fwi%forward,fwi%nb_elem,fwi%DoF,fwi%time_scheme, &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%k_max,           &
         fwi%epsilon,fwi%source_loc,fwi%receiver_loc)

    fwi%n_time_step=fwi%forward%n_time_step
    print*,'Il y a ',fwi%n_time_step,'time steps'
    fwi%dt=fwi%forward%dt
    fwi%dx=fwi%forward%dx
    call Sparse2Full(fwi%forward%M,fwi%M)
    call Sparse2Full(fwi%forward%Ap,fwi%Ap)
    call Sparse2Full(fwi%forward%App,fwi%App)
    call Sparse2Full(fwi%forward%Av,fwi%Av)
    
    allocate(fwi%P(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%U(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%QP(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%QU(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%P_received(0:fwi%n_time_step,2))
    allocate(fwi%U_received(0:fwi%n_time_step,2))
    fwi%P(:,0)=0
    fwi%QP(:,0)=0
    fwi%U(:,0)=0
    fwi%QU(:,0)=0
    fwi%P_received(:,:)=0.0
    fwi%U_received(:,:)=0.0
    
    if (fwi%adjoint_test) then
       ! allocate(fwi%P(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
       ! allocate(fwi%U(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
       allocate(fwi%FP(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
       allocate(fwi%FU(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
       allocate(fwi%FP_half(fwi%DoF*fwi%nb_elem,1:fwi%n_time_step))
       allocate(fwi%FU_half(fwi%DoF*fwi%nb_elem,1:fwi%n_time_step))
       allocate(fwi%GP(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
       allocate(fwi%GU(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    end if

    current_time_step=0
    t=0
    do current_time_step=1,fwi%n_time_step
       t=current_time_step*fwi%forward%dt
       call one_forward_step(fwi%forward,t)
       fwi%P(:,current_time_step)=fwi%forward%P
       fwi%U(:,current_time_step)=fwi%forward%U
       fwi%P_received(current_time_step,1)=t
       fwi%P_received(current_time_step,2)=fwi%forward%P(fwi%DoF*(fwi%receiver_loc-1)+1)
       fwi%U_received(current_time_step,1)=t
       fwi%U_received(current_time_step,2)=fwi%forward%U(fwi%DoF*(fwi%receiver_loc-1)+1)
       
       if (fwi%adjoint_test) then
          fwi%FP(:,current_time_step)=fwi%forward%RHSp
          fwi%FU(:,current_time_step)=fwi%forward%RHSv
          if (fwi%time_scheme.eq.'RK4') then
             fwi%FP_half(:,current_time_step)=fwi%forward%RHSp_half
             fwi%FU_half(:,current_time_step)=fwi%forward%RHSv_half
          end if
       end if
    end do

    if (fwi%adjoint_test) then
       fwi%FP(:,fwi%forward%n_time_step)=0.0
       fwi%FU(:,fwi%forward%n_time_step)=0.0
       fwi%FP(:,0)=0.0
       fwi%FU(:,0)=0.0
       call fwi_extract_g(fwi)
    end if

    if (fwi%animation) then
       do i=0,fwi%n_time_step
          if (modulo(i,10).eq.0) then
             call print_vect(fwi%P(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'P')
          else if (i.eq.0) then
             call print_vect(fwi%P(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'P')
          end if
       end do
    end if

    if (fwi%adjoint_test) then
       test%M=fwi%forward%M
       test%Minv=fwi%forward%Minv
       test%Ap=fwi%forward%Ap
       test%Av=fwi%forward%Av
       test%App=fwi%forward%App
    end if
    
    call free_acoustic_problem(fwi%forward)
    
    ! print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    ! print*,'%%%%%%%%%%%%%%%%%%%%% backward %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    ! print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call init_adjoint_problem(fwi%backward,fwi%nb_elem,fwi%DoF,fwi%time_scheme, &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%k_max,           &
         fwi%epsilon,fwi%source_loc,fwi%receiver_loc,fwi%data_P,fwi%data_U,     &
         fwi%P_received,fwi%U_received,fwi%strategy,fwi%scalar_product)

   ! print*,'test n_time_step'
    if (fwi%forward%n_time_step.ne.fwi%backward%n_time_step) then
       print*,'ERROR n_time_step'
       print*,'forward : ',fwi%forward%n_time_step,fwi%forward%dt
       print*,'backward : ',fwi%backward%n_time_step,fwi%backward%dt
       STOP
    end if

    fwi%n_time_step=fwi%backward%n_time_step

    if (fwi%adjoint_test) then
       allocate(fwi%DP(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
       allocate(fwi%DU(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    end if
    
    fwi%QP(:,fwi%n_time_step)=0
    fwi%QU(:,fwi%n_time_step)=0

    current_time_step=fwi%n_time_step
    t=0

    do current_time_step=fwi%n_time_step-1,0,-1
       !do current_time_step=1,fwi%n_time_step
       t=current_time_step*fwi%backward%dt
       call one_backward_step(fwi%backward,t)
       fwi%QP(:,current_time_step)=fwi%backward%P
       fwi%QU(:,current_time_step)=fwi%backward%U


       if (fwi%adjoint_test) then
          if (fwi%time_scheme.eq.'RK4') then
             fwi%DP(:,current_time_step)=fwi%backward%RHSp
             fwi%DU(:,current_time_step)=fwi%backward%RHSv
          else if (fwi%time_scheme.eq.'AB3') then
             fwi%DP(:,current_time_step)=fwi%backward%RHSp
             fwi%DU(:,current_time_step)=fwi%backward%RHSv
          end if
       end if
    end do

    if (fwi%adjoint_test) then
       fwi%DP(:,0)=0.0
       fwi%DU(:,0)=0.0
       fwi%DP(:,fwi%backward%n_time_step)=0.0
       fwi%DU(:,fwi%backward%n_time_step)=0.0
    end if

    
    if (fwi%animation) then
       ! print*,'Calculus done'
       do i=0,fwi%n_time_step
          if (modulo(i,10).eq.0) then
             call print_vect(fwi%QP(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'QP')
          else if (i.eq.0) then
             call print_vect(fwi%QP(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'QP')
          end if
       end do
    end  if

    
    ! RTM=0.0
    ! do i=1,size(RTM)
    !    do j=0,fwi%n_time_step
    !       RTM(i)=RTM(i)+fwi%P(i,j)*fwi%QP(i,fwi%n_time_step-j)
    !    end do
    ! end do

    ! open(unit=23,file='RTM.dat')
    ! do i=1,size(RTM)
    !    write(23,*) RTM(i)
    ! end do
    ! close(23)
    
    if (fwi%adjoint_test) then
       print*,'%%%%%%% ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
       print*,'Inner Product PU/DPDU',inner_product(fwi%P,fwi%U,fwi%DP,fwi%DU)
       print*,'Inner Product QPQU/GPGU',inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU)
       print*,'Difference relative :',(inner_product(fwi%P,fwi%U,fwi%DP,fwi%DU)-   &
            inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU))/&
            inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU)

       print*,'test',size(fwi%P,1),size(fwi%P,2)

       print*,'M Inner Product PU/DPDU',inner_product_M(fwi%P,fwi%U,fwi%DP,fwi%DU,test%M)
       print*,'M Inner Product QPQU/GPGU',inner_product_M(fwi%QP,fwi%QU,fwi%GP,fwi%GU,test%M)
       print*,'Difference relative :',(inner_product_M(fwi%P,fwi%U,fwi%DP,fwi%DU,test%M)-   &
                                       inner_product_M(fwi%QP,fwi%QU,fwi%GP,fwi%GU,test%M))/&
                                       inner_product_M(fwi%QP,fwi%QU,fwi%GP,fwi%GU,test%M)
       print*,'%%%%%% END ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    end if
    
    call free_adjoint_problem(fwi%backward)
    
    ! print*,'test gradient'
    call CPU_TIME(t1)    
    if (fwi%strategy.eq.'ATD') then
       call gradient_ATD_computation(fwi)
    else if (fwi%strategy.eq.'DTA') then
       call gradient_DTA_computation(fwi)
    end if
    call CPU_TIME(t2)
    print*,'Temps de calcul du gradient est de : ',t2-t1

    call update_model(fwi)

    
    print*,'%%%% COST FUNCTION %%%%%'
    call cost_function(fwi%P_received,fwi%data_P,fwi%final_time,fwi%current_iter)
    deallocate(fwi%P)
    deallocate(fwi%U)
    deallocate(fwi%QP)
    deallocate(fwi%QU)
    deallocate(fwi%P_received)
    deallocate(fwi%U_received)
    deallocate(fwi%M)
    deallocate(fwi%Ap)
    deallocate(fwi%App)
    deallocate(fwi%Av)
    if (fwi%adjoint_test) then
       deallocate(fwi%FP)
       deallocate(fwi%FU)
       deallocate(fwi%FP_half)
       deallocate(fwi%FU_half)
       deallocate(fwi%GP)
       deallocate(fwi%GU)
    end if

    
     ! print*,'%%%%%%%%%%%%%%%%%% ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%'
     ! call  init_adjoint_test(test,fwi%n_time_step,fwi%time_scheme,fwi%backward%dt,&
     !      test%Ap,test%Av,test%App,test%M,test%Minv,&
     !      fwi%nb_elem,fwi%DoF,fwi%backward%dx,fwi%FP,fwi%FU,fwi%FP_half,fwi%FU_half,fwi%DP,fwi%DU)
     ! call forward_test(test)
     ! call extract_g(test)
     ! call backward_test(test)
     ! print*,'Inner Product PU/DPDU',inner_product(test%P,test%U,test%DP,test%DU)
     ! print*,'Inner Product QPQU/GPGU',inner_product(test%QP,test%QU,test%GP,test%GU)
     ! print*,'M Inner Product PU/DPDU',inner_product_M(test%P,test%U,test%DP,test%DU,test%M)
     ! print*,'M Inner Product QPQU/GPGU',inner_product_M(test%QP,test%QU,test%GP,test%GU,test%M)
     ! print*,'norm test P',(norm_test(fwi%P)-norm_test(test%P))/norm_test(test%P)
     ! print*,'norm test U',(norm_test(fwi%U)-norm_test(test%U))/norm_test(test%U)
     ! print*,'norm test QP',(norm_test(fwi%QP)-norm_test(test%QP))/norm_test(test%QP)
     ! print*,'norm test QU',(norm_test(fwi%QU)-norm_test(test%QU))/norm_test(test%QU)
     ! print*,'norm test FP',norm_test(fwi%FP)-norm_test(test%FP)
     ! print*,'norm test FU',norm_test(fwi%FU)-norm_test(test%FU)
     ! print*,'norm test FP_half',norm_test(fwi%FP_half)-norm_test(test%FP_half)
     ! print*,'norm test FU_half',norm_test(fwi%FU_half)-norm_test(test%FU_half)
     ! print*,'norm test DP',norm_test(fwi%DP)-norm_test(test%DP)     
     ! print*,'norm test DU',norm_test(fwi%DU)-norm_test(test%DU)   
     
     ! print*,'%%%%%%%%%%%%%%%% END ADJOINT TEST %%%%%%%%%%%%%%%%%'

    !--------------------- Animation ----------------------------------------------
  !n_frame=int(forward%n_time_step/forward%n_display)
     ! open(unit=78,file='script.gnuplot',action='write')
     ! write(78,*)'load "trace1.gnuplot"'
     ! write(78,*)'n=',fwi%n_time_step
     ! write(78,*)'a=',30
     ! write(78,*)'load "trace2.gnuplot"'
     ! close(78)
  
     ! call system('gnuplot script.gnuplot')
    
    end subroutine one_fwi_step
  
  subroutine fwi_extract_g(fwi)
    type(t_fwi),intent(inout) :: fwi

    real,dimension(size(fwi%P,1)) :: Pk1,Pk2
    real,dimension(size(fwi%P,1)) :: Uk1,Uk2
    integer :: i
    real,dimension(size(fwi%P,1)) :: zero
    zero=0.0
    fwi%GP=0.0
    fwi%GU=0.0

    if (fwi%time_scheme.eq.'LF') then
       print*,'not possible yet for Leap Frog'
       STOP


    else if (fwi%time_scheme.eq.'RK4') then
       do i=1,fwi%n_time_step
          call RK4_forward(fwi%GP(:,i),fwi%GU(:,i),                             &
               fwi%forward%Ap,fwi%forward%Av,fwi%forward%App,                   &
               fwi%FP(:,i-1),fwi%FP_half(:,i),fwi%FP(:,i),                      &
               fwi%FU(:,i-1),fwi%FU_half(:,i),fwi%FU(:,i))
       end do

    else if (fwi%time_scheme.eq.'AB3') then

       Pk1=0.0
       Pk2=0.0
       Uk1=0.0
       Uk2=0.0

       call  AB3_forward(fwi%GP(:,0),fwi%GU(:,0),fwi%forward%Ap,fwi%forward%Av, &
                         fwi%forward%App,fwi%FP(:,0),fwi%FU(:,0),               &
                         zero,zero,                                             &
                         zero,zero,                                             &
                         Pk1,Pk2,Uk1,Uk2)
       
       
       call  AB3_forward(fwi%GP(:,1),fwi%GU(:,1),fwi%forward%Ap,fwi%forward%Av, &
                         fwi%forward%App,fwi%FP(:,1),fwi%FU(:,1),               &
                         fwi%FP(:,0),fwi%FU(:,0),                               &
                         zero,zero,                                             &
                         Pk1,Pk2,Uk1,Uk2)


       
       do i=2,fwi%n_time_step
          call  AB3_forward(fwi%GP(:,i),fwi%GU(:,i),fwi%forward%Ap,             &
                            fwi%forward%Av,fwi%forward%App,                     &
                            fwi%FP(:,i),fwi%FU(:,i),                            &
                            fwi%FP(:,i-1),fwi%FU(:,i-1),                        &
                            fwi%FP(:,i-2),fwi%FU(:,i-2),                        &
                            Pk1,Pk2,Uk1,Uk2)
       end do

       
    else
       print*,'not recongnized time scheme'
    end if

  end subroutine fwi_extract_g



  subroutine gradient_ATD_computation(fwi)
    type(t_fwi),intent(inout) :: fwi
    real,dimension(fwi%DoF)       :: dtu,dtp,qu,qp
    real,dimension(fwi%DoF*fwi%nb_elem) :: u0,u1,u2
    real,dimension(fwi%DoF*fwi%nb_elem) :: p0,p1,p2
    real,dimension(fwi%DoF*fwi%nb_elem) :: dtp_glob

    integer :: nb_elem,DoF,beg_node,end_node
    integer :: i,j,k
    real    :: moyenne


    
    nb_elem=fwi%nb_elem   ! For sake of
    DoF=fwi%DoF           ! readibility
    

    fwi%gradient_vp=0.0
    fwi%gradient_rho=0.0
    
    
    do j=0,fwi%n_time_step
       do i=1,nb_elem
          if (j.eq.0) then
             beg_node=(i-1)*DoF+1
             end_node=i*DoF
             dtu=(fwi%U(beg_node:end_node,j+1)-fwi%U(beg_node:end_node,j))/(fwi%dt)
             dtp=(fwi%P(beg_node:end_node,j+1)-fwi%P(beg_node:end_node,j))/(fwi%dt)
             qp=fwi%QP(beg_node:end_node,j)
             qu=fwi%QU(beg_node:end_node,j)

          else if (j.eq.fwi%n_time_step) then
             beg_node=(i-1)*DoF+1
             end_node=i*DoF
             dtu=(fwi%U(beg_node:end_node,j)-fwi%U(beg_node:end_node,j-1))/(fwi%dt)
             dtp=(fwi%P(beg_node:end_node,j)-fwi%P(beg_node:end_node,j-1))/(fwi%dt)
             qp=fwi%QP(beg_node:end_node,j)
             qu=fwi%QU(beg_node:end_node,j)
          else
             beg_node=(i-1)*DoF+1
             end_node=i*DoF
             dtu=(fwi%U(beg_node:end_node,j+1)-fwi%U(beg_node:end_node,j-1))/(2*fwi%dt)
             dtp=(fwi%P(beg_node:end_node,j+1)-fwi%P(beg_node:end_node,j-1))/(2*fwi%dt)
             qp=fwi%QP(beg_node:end_node,j)
             qu=fwi%QU(beg_node:end_node,j)
          end if


          if (fwi%bernstein) then

             dtp=matmul(B2L,dtp)
             dtu=matmul(B2L,dtu)
             qp=matmul(B2L,qp)
             qu=matmul(B2L,qu)

          end if

          moyenne=0.0
          do k=1,DoF
             moyenne=dtp(k)*qp(k)
          end do
          moyenne=moyenne/DoF

          fwi%gradient_vp(i)=fwi%gradient_vp(i)                              &
               -2.0/(fwi%density_model(i)*fwi%velocity_model(i)**3)*moyenne

       end do
    end do

       fwi%gradient_vp=fwi%gradient_vp*fwi%dx*fwi%dt        

     end subroutine gradient_ATD_computation

     subroutine gradient_DTA_computation(fwi)
       type(t_fwi),intent(inout) :: fwi
       real,dimension(fwi%DoF)       :: dtu,dtp,qu,qp
       real,dimension(fwi%DoF*fwi%nb_elem) :: u0,u1,u2
       real,dimension(fwi%DoF*fwi%nb_elem) :: p0,p1,p2
       real,dimension(fwi%DoF*fwi%nb_elem) :: dtp_glob

       integer :: nb_elem,DoF,beg_node,end_node
       integer :: i,j,k
       real    :: one_parameter,moyenne
       real,dimension(fwi%DoF,fwi%DoF) :: M_loc
       real,dimension(fwi%DoF,fwi%DoF*fwi%nb_elem) :: Ap_loc
       real,dimension(fwi%DoF,fwi%DoF*fwi%nb_elem) :: App_loc
       real,dimension(fwi%DoF,fwi%DoF*fwi%nb_elem) :: Av_loc

       type(sparse_matrix)                         :: Ap_sparse_loc
       type(sparse_matrix)                         :: App_sparse_loc
       type(sparse_matrix)                         :: Av_sparse_loc


       nb_elem=fwi%nb_elem   ! For sake of
       DoF=fwi%DoF           ! readibility
    
       
       fwi%gradient_vp=0.0
       fwi%gradient_rho=0.0


       if (fwi%time_scheme.eq.'RK4') then
          do j=0,fwi%n_time_step
             do i=1,nb_elem
                M_loc=fwi%M((i-1)*DoF+1:i*DoF,(i-1)*DoF+1:i*DoF)

             if (j.eq.0) then
                beg_node=(i-1)*DoF+1
                end_node=i*DoF
                dtu=(fwi%U(beg_node:end_node,j+1)-fwi%U(beg_node:end_node,j))/(fwi%dt)
                dtp=(fwi%P(beg_node:end_node,j+1)-fwi%P(beg_node:end_node,j))/(fwi%dt)
                qp=fwi%QP(beg_node:end_node,j)
                qu=fwi%QU(beg_node:end_node,j)

             else if (j.eq.fwi%n_time_step) then
                beg_node=(i-1)*DoF+1
                end_node=i*DoF
                dtu=(fwi%U(beg_node:end_node,j)-fwi%U(beg_node:end_node,j-1))/(fwi%dt)
                dtp=(fwi%P(beg_node:end_node,j)-fwi%P(beg_node:end_node,j-1))/(fwi%dt)
                qp=fwi%QP(beg_node:end_node,j)
                qu=fwi%QU(beg_node:end_node,j)
             else
                beg_node=(i-1)*DoF+1
                end_node=i*DoF
                dtu=(fwi%U(beg_node:end_node,j+1)-fwi%U(beg_node:end_node,j-1))/(2*fwi%dt)
                dtp=(fwi%P(beg_node:end_node,j+1)-fwi%P(beg_node:end_node,j-1))/(2*fwi%dt)
                qp=fwi%QP(beg_node:end_node,j)
                qu=fwi%QU(beg_node:end_node,j)
             end if

             if (fwi%strategy.eq.'DTA') then
                dtp=matmul(M_loc,dtp)
                dtu=matmul(M_loc,dtu)
                if (fwi%scalar_product.eq.'M') then
                   dtp=matmul(M_loc,dtp)
                   dtu=matmul(M_loc,dtu)
                end if
             end if

             moyenne=0.0
             do k=1,DoF
                moyenne=dtp(k)*qp(k)
             end do
             moyenne=moyenne/DoF

             fwi%gradient_vp(i)=fwi%gradient_vp(i)                              &
                  -2.0/(fwi%density_model(i)*fwi%velocity_model(i)**3)*moyenne
          end do
       end do

    else

       fwi%gradient_vp=fwi%gradient_vp*fwi%dx*fwi%dt        

       if ((fwi%strategy.eq.'DTA').and.(fwi%time_scheme.eq.'AB3')) then

       print*,'je suis passe par là'

          fwi%gradient_vp=0.0

          fwi%App(1,1)=fwi%App(1,1)/fwi%velocity_model(1)
          fwi%App(fwi%DoF*fwi%nb_elem,fwi%DoF*fwi%nb_elem)=                     &
               fwi%App(fwi%DoF*fwi%nb_elem,fwi%DoF*fwi%nb_elem)/                &
               fwi%velocity_model(fwi%DoF)
          

             do i=1,nb_elem
                Ap_loc=fwi%Ap((i-1)*DoF+1:i*DoF,:)/(0.5*fwi%velocity_model(i))
                App_loc=fwi%Ap((i-1)*DoF+1:i*DoF,:)
                Av_loc=fwi%Ap((i-1)*DoF+1:i*DoF,:)
                
                beg_node=(i-1)*DoF+1
                end_node=i*DoF

                call Full2Sparse(Ap_loc,Ap_sparse_loc)
                call Full2Sparse(App_loc,App_sparse_loc)
                call Full2Sparse(Av_loc,Av_sparse_loc)

                do j=0,fwi%n_time_step
                qp=fwi%QP(beg_node:end_node,j)


                u0=(fwi%U(:,j))
                p0=(fwi%P(:,j))
                if (j.eq.0) then
                   u1=0.0
                   u2=0.0
                   p1=0.0
                   p2=0.0
                else if (j.eq.1) then
                   u1=(fwi%U(:,j-1))
                   p1=(fwi%P(:,j-1))
                   u2=0.0
                   p2=0.0
                else
                   u1=(fwi%U(:,j-1))
                   p1=(fwi%P(:,j-1))
                   u2=(fwi%U(:,j-2))
                   p2=(fwi%P(:,j-2))

                   dtp=-23.0/12.0*(sparse_matmul(App_sparse_loc,p0)+sparse_matmul(Ap_sparse_loc,u0))  &
                        +16.0/12.0*(sparse_matmul(App_sparse_loc,p1)+sparse_matmul(Ap_sparse_loc,u1)) &
                        -5.0/12.0*(sparse_matmul(App_sparse_loc,p2)+sparse_matmul(Ap_sparse_loc,u2))

                   moyenne=0.0
                   do k=1,fwi%DoF
                      moyenne=dtp(k)*qp(k)
                   end do

                   fwi%gradient_vp(i)=fwi%gradient_vp(i)+moyenne
                end if
             end do
          end do
       end if
    end if

       
     end subroutine gradient_DTA_computation


     subroutine update_model(fwi)
       type(t_fwi),intent(inout) :: fwi

       integer                         :: i,j
       integer,parameter               :: nb_section=100
       integer,dimension(nb_section,2) :: lim_section
       real,dimension(nb_section)      :: grad
       real,parameter                  :: beta=0.05
       character(len=50)               :: fichier


       grad=0.0

       do i=1,nb_section
          lim_section(i,1)=(i-1)*int(fwi%nb_elem/nb_section)+1
          lim_section(i,2)=i*int(fwi%nb_elem/nb_section)
       end do
       lim_section(nb_section,2)=fwi%nb_elem

       !grad(1:33)=0.0
       do i=1,nb_section
          grad(i)=sum(fwi%gradient_vp(lim_section(i,1):lim_section(i,2)))
       end do
       
       if (maxval(abs(grad)).eq.0.0) then
          print*,'Le gradient est nul, fin de l algo'
          STOP
       else
          !grad=(1.0/maxval(abs(grad)))*grad
          grad=(1.0/grad(50))*grad
       end if

       open(unit=19,file='grad.dat')
       do i=1,size(grad)
          write(19,*) i,grad(i)
       end do

       do i=1,nb_section
          do j=lim_section(i,1),lim_section(i,2)
             fwi%velocity_model(j)=fwi%velocity_model(j)-beta*grad(i)
             if (fwi%velocity_model(j).lt.1.0) then
                fwi%velocity_model(j)=1.0
             else if (fwi%velocity_model(j).gt.2.0) then
                fwi%velocity_model(j)=2.0
             end if
          end do
       end do

       write(fichier,"(A,A,I0,'.dat')") "../Bernstein/Files/",'VP',fwi%current_iter
       open(unit=28,file=fichier)
       do i=1,fwi%nb_elem
          write(28,*) i,fwi%velocity_model(i)
       end do

       ! do i=1,nb_section
       !    print*,'model test',fwi%velocity_model(lim_section(i,1))
       ! end do
       
     end subroutine update_model


  subroutine cost_function(P_received,data_P,final_time,iter)
    real,dimension(:,:),intent(in) :: P_received
    real,dimension(:,:),intent(in) :: data_P
    real               ,intent(in) :: final_time
    integer            ,intent(in) :: iter
    real                           :: CF
    real :: eval_P_received,eval_data_P
    real :: dt,t
    integer :: j,i,n_time_step

    dt=P_received(2,1)-P_received(1,1)
    dt=min(dt,data_P(2,1)-data_P(1,1))

    n_time_step=max(size(P_received,1),size(data_P,1))

    CF=0.0
    
    do j=2,n_time_step-1
       t=(j-1)*dt
       eval_P_received=0.0
       eval_data_P=0.0

       i=1
       do while ((P_received(i,1).le.t))
          i=i+1
       end do
       i=i-1

       if (i.eq.0) then
          print*,'test',i,t,P_received(1,1)
       end if
       
       eval_P_received=P_received(i+1,2)*(t-P_received(i,1))                    &
                      +P_received(i,2)*(P_received(i+1,1)-t)
       eval_P_received=eval_P_received/(P_received(i+1,1)-P_received(i,1))


       i=1
       do while ((data_P(i,1).le.t))
          i=i+1
       end do
       i=i-1

       
       if (i.eq.0) then
          print*,'test',i,t,data_P(1,1)
       end if
       

       eval_data_P=data_P(i+1,2)*(t-data_P(i,1))                    &
            +data_P(i,2)*(data_P(i+1,1)-t)
       eval_data_P=eval_data_P/(data_P(i+1,1)-data_P(i,1))

       CF=CF+(eval_P_received-eval_data_P)**2

       write(26,*) t,eval_data_P
       write(27,*) t,eval_P_received
       
    end do
    CF=CF*dt
    CF=0.5*CF

    print*,'CF=',CF
    write(33,*) iter,CF

    
    if (CF.eq.0.0) then
       print*,'Le gradient a parfaitement convergé'
       STOP
    end if

       
  end subroutine cost_function


  

end module m_fwi
