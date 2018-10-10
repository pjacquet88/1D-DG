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
     integer                         :: model_size
     real,dimension(:)  ,allocatable :: gradJ_velocity
     real,dimension(:)  ,allocatable :: gradJ_density
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
     
     integer                         :: nb_elem
     integer                         :: DoF
     character(len=20)               :: time_scheme
     real                            :: total_length
     real                            :: final_time
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
  end type t_fwi
    

contains

  subroutine init_fwi(fwi,nb_iter,velocity_model,density_model,data_P,data_U,   &
                      nb_elem,DoF,time_scheme,total_length,final_time,alpha,    &
                      bernstein,k_max,epsilon,source_loc,receiver_loc,strategy, &
                      scalar_product)
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
    
    fwi%nb_iter=nb_iter
    fwi%model_size=size(density_model)

    allocate(fwi%velocity_model(fwi%model_size))
    allocate(fwi%density_model(fwi%model_size))

    fwi%velocity_model=velocity_model
    fwi%density_model=density_model

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
    
    allocate(fwi%gradJ_velocity(fwi%model_size))
    allocate(fwi%gradJ_density(fwi%model_size))

    fwi%signal='flat'
    fwi%boundaries='ABC'
  end subroutine init_fwi
  

  subroutine free_fwi(fwi)
    type(t_fwi),intent(inout) :: fwi

    deallocate(fwi%velocity_model)
    deallocate(fwi%density_model)
    deallocate(fwi%data_P)
    deallocate(fwi%data_U)
    deallocate(fwi%gradJ_velocity)
    deallocate(fwi%gradJ_density)
  end subroutine free_fwi


  subroutine one_fwi_step(fwi)
    type(t_fwi),intent(inout) :: fwi
    integer                   :: current_time_step
    real                      :: t

    real,dimension(fwi%DoF*fwi%nb_elem) :: RTM
    integer                             :: i,j
    type(t_adjoint_test)                :: test

    !------------------------------ Forward -------------------------------------
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%% FORWARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call init_acoustic_problem(fwi%forward,fwi%nb_elem,fwi%DoF,fwi%time_scheme, &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%k_max,           &
         fwi%epsilon,fwi%source_loc,fwi%receiver_loc)

    fwi%n_time_step=fwi%forward%n_time_step
 
    allocate(fwi%P(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%U(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%FP(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%FU(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%FP_half(fwi%DoF*fwi%nb_elem,1:fwi%n_time_step))
    allocate(fwi%FU_half(fwi%DoF*fwi%nb_elem,1:fwi%n_time_step))
    allocate(fwi%GP(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%GU(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    
    allocate(fwi%P_received(0:fwi%n_time_step,2))
    allocate(fwi%U_received(0:fwi%n_time_step,2))
    fwi%P(:,0)=0
    fwi%U(:,0)=0

    current_time_step=0
    t=0
    do current_time_step=1,fwi%n_time_step
       t=current_time_step*fwi%forward%dt
       call one_forward_step(fwi%forward,t)
       fwi%P(:,current_time_step)=fwi%forward%P
       fwi%U(:,current_time_step)=fwi%forward%U
       fwi%P_received(current_time_step,1)=t
       fwi%P_received(current_time_step,2)=fwi%forward%P(fwi%receiver_loc)
       fwi%U_received(current_time_step,1)=t
       fwi%U_received(current_time_step,2)=fwi%forward%U(fwi%receiver_loc)
       fwi%FP(:,current_time_step)=fwi%forward%RHSp
       fwi%FU(:,current_time_step)=fwi%forward%RHSv
       if (fwi%time_scheme.eq.'RK4') then
       fwi%FP_half(:,current_time_step)=fwi%forward%RHSp_half
       fwi%FU_half(:,current_time_step)=fwi%forward%RHSv_half
       end if
    end do
    fwi%FP(:,fwi%forward%n_time_step)=0.0
    fwi%FU(:,fwi%forward%n_time_step)=0.0
    ! fwi%FP(:,0)=0.0
    ! fwi%FU(:,0)=0.0

   call fwi_extract_g(fwi)

    print*,'Calculus done'
    do i=0,fwi%n_time_step
       if (modulo(i,20).eq.0) then
          call print_vect(fwi%P(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'P')
       else if (i.eq.0) then
          call print_vect(fwi%P(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'P')
       end if
    end do

    test%M=fwi%forward%M
    test%Minv=fwi%forward%Minv
    test%Ap=fwi%forward%Ap
    test%Av=fwi%forward%Av
    test%App=fwi%forward%App

    call free_acoustic_problem(fwi%forward)
    print*,'Initialization done'
    
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%% backward %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call init_adjoint_problem(fwi%backward,fwi%nb_elem,fwi%DoF,fwi%time_scheme, &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%k_max,           &
         fwi%epsilon,fwi%source_loc,fwi%receiver_loc,fwi%data_P,fwi%data_U,     &
         fwi%P_received,fwi%U_received,fwi%strategy,fwi%scalar_product)

    print*,'test n_time_step'
    if (fwi%forward%n_time_step.ne.fwi%backward%n_time_step) then
       print*,'ERROR n_time_step'
       print*,'forward : ',fwi%forward%n_time_step,fwi%forward%dt
       print*,'backward : ',fwi%backward%n_time_step,fwi%backward%dt
       STOP
    end if

    fwi%n_time_step=fwi%backward%n_time_step

    print*,'ALLOCATION'
    allocate(fwi%QP(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%QU(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%DP(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%DU(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    fwi%QP(:,fwi%n_time_step)=0
    fwi%QU(:,fwi%n_time_step)=0
    print*,'FIN ALLOCATION'
    current_time_step=fwi%n_time_step
    t=0
    print*,'DEBUT COPIE'
    do current_time_step=fwi%n_time_step-1,0,-1
       !do current_time_step=1,fwi%n_time_step
       t=current_time_step*fwi%backward%dt
       call one_backward_step(fwi%backward,t)
       fwi%QP(:,current_time_step)=fwi%backward%P
       fwi%QU(:,current_time_step)=fwi%backward%U
       if (fwi%time_scheme.eq.'RK4') then
          fwi%DP(:,current_time_step)=fwi%backward%RHSp
          fwi%DU(:,current_time_step)=fwi%backward%RHSv
       else if (fwi%time_scheme.eq.'AB3') then
          fwi%DP(:,current_time_step)=fwi%backward%RHSp
          fwi%DU(:,current_time_step)=fwi%backward%RHSv
       end if
    end do
!        fwi%DP(:,0)=0.0
!        fwi%DU(:,0)=0.0
    fwi%DP(:,fwi%backward%n_time_step)=0.0
    fwi%DU(:,fwi%backward%n_time_step)=0.0
    print*,'FIN COPIE'

    print*,'Calculus done'
    do i=0,fwi%n_time_step
       if (modulo(i,20).eq.0) then
          call print_vect(fwi%QP(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'QP')
       else if (i.eq.0) then
          call print_vect(fwi%QP(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'QP')
       end if
    end do

    RTM=0.0
    do i=1,size(RTM)
       do j=0,fwi%n_time_step
          RTM(i)=RTM(i)+fwi%P(i,j)*fwi%QP(i,fwi%n_time_step-j)
       end do
    end do

    open(unit=23,file='RTM.dat')
    do i=1,size(RTM)
       write(23,*) RTM(i)
    end do
    close(23)


    print*,'%%%%%%% ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'Inner Product PU/DPDU',inner_product(fwi%P,fwi%U,fwi%DP,fwi%DU)
    print*,'Inner Product QPQU/GPGU',inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU)
    print*,'Difference relative :',(inner_product(fwi%P,fwi%U,fwi%DP,fwi%DU)-   &
                                    inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU))/&
                                    inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU)

        print*,'M Inner Product PU/DPDU',inner_product_M(fwi%P,fwi%U,fwi%DP,fwi%DU,fwi%backward%M)
    print*,'M Inner Product QPQU/GPGU',inner_product_M(fwi%QP,fwi%QU,fwi%GP,fwi%GU,fwi%backward%M)
    print*,'Difference relative :',(inner_product_M(fwi%P,fwi%U,fwi%DP,fwi%DU,fwi%backward%M)-   &
                                    inner_product_M(fwi%QP,fwi%QU,fwi%GP,fwi%GU,fwi%backward%M))/&
                                    inner_product_M(fwi%QP,fwi%QU,fwi%GP,fwi%GU,fwi%backward%M)

    print*,'test',maxval(fwi%backward%M%Values),minval(fwi%backward%M%Values)
    print*,'test',maxval(test%M%Values),minval(test%M%Values)
    print*,'%%%%%% END ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'


     print*,'%%%%%%%%%%%%%%%%%% ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%'
     call  init_adjoint_test(test,fwi%n_time_step,fwi%time_scheme,fwi%backward%dt,&
          test%Ap,test%Av,test%App,test%M,test%Minv,&
          fwi%nb_elem,fwi%DoF,fwi%backward%dx,fwi%FP,fwi%FU,fwi%FP_half,fwi%FU_half,fwi%DP,fwi%DU)
     call forward_test(test)
     call extract_g(test)
     call backward_test(test)
     print*,'Inner Product PU/DPDU',inner_product(test%P,test%U,test%DP,test%DU)
     print*,'Inner Product QPQU/GPGU',inner_product(test%QP,test%QU,test%GP,test%GU)
     print*,'M Inner Product PU/DPDU',inner_product_M(test%P,test%U,test%DP,test%DU,test%M)
     print*,'M Inner Product QPQU/GPGU',inner_product_M(test%QP,test%QU,test%GP,test%GU,test%M)
     print*,'norm test P',(norm_test(fwi%P)-norm_test(test%P))/norm_test(test%P)
     print*,'norm test U',(norm_test(fwi%U)-norm_test(test%U))/norm_test(test%U)
     print*,'norm test QP',(norm_test(fwi%QP)-norm_test(test%QP))/norm_test(test%QP)
     print*,'norm test QU',(norm_test(fwi%QU)-norm_test(test%QU))/norm_test(test%QU)
     print*,'norm test FP',norm_test(fwi%FP)-norm_test(test%FP)
     print*,'norm test FU',norm_test(fwi%FU)-norm_test(test%FU)
     print*,'norm test FP_half',norm_test(fwi%FP_half)-norm_test(test%FP_half)
     print*,'norm test FU_half',norm_test(fwi%FU_half)-norm_test(test%FU_half)
     print*,'norm test DP',norm_test(fwi%DP)-norm_test(test%DP)     
     print*,'norm test DU',norm_test(fwi%DU)-norm_test(test%DU)     
     print*,'%%%%%%%%%%%%%%%% END ADJOINT TEST %%%%%%%%%%%%%%%%%'
     
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

end module m_fwi
