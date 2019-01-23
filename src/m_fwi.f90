module m_fwi
  use m_kind
  use m_file_function
  use m_matrix
  use m_time_scheme
  use m_acoustic
  use m_adjoint
  use m_adjoint_test
  implicit none

  type t_fwi
     type(acoustic_problem)              :: forward         ! the forward state
     type(adjoint_problem)               :: backward        ! the backwardstate 
     integer                             :: nb_iter         ! nb iteration of the fwi
     integer                             :: current_iter    ! nb of the current iteration
     real(mp),dimension(:)  ,allocatable :: velocity_model  ! the velocity model
     real(mp),dimension(:)  ,allocatable :: density_model   ! the density model
     integer                             :: n_time_step     ! nb of time step
     real(mp),dimension(:,:),allocatable :: data_P          ! pressure measured at the receiver
     real(mp),dimension(:,:),allocatable :: data_U          ! velocity measured at the reciever
     real(mp),dimension(:,:),allocatable :: P_received      ! pressure simulated at the receiver
     real(mp),dimension(:,:),allocatable :: U_received      ! velocity simulated at the receiver
     real(mp),dimension(:,:),allocatable :: P               ! forward pressure state in space and time
     real(mp),dimension(:,:),allocatable :: U               ! forward velocity state in space and time
     real(mp),dimension(:,:),allocatable :: QP              ! backward pressure state in space and time
     real(mp),dimension(:,:),allocatable :: QU              ! backward velocity state in space and time
     real(mp),dimension(:,:),allocatable :: FP              ! forward dtP equation RHS in space and time
     real(mp),dimension(:,:),allocatable :: FU              ! forward dtU equation RHS in space and time
     real(mp),dimension(:,:),allocatable :: FP_half         ! forward dtP equation in space and time (on half time)
     real(mp),dimension(:,:),allocatable :: FU_half         ! forward dtU equation in space and time (on half time)
     real(mp),dimension(:,:),allocatable :: DP              ! backward dtQP equation in space and time
     real(mp),dimension(:,:),allocatable :: DU              ! backward dtQU equation in space and time
     real(mp),dimension(:,:),allocatable :: GP              ! combined forward dtP RHS in space and time (G=EF)
     real(mp),dimension(:,:),allocatable :: GU              ! combined foward dtU RHS in space and time (G=EF)
     real(mp),dimension(:)  ,allocatable :: gradient_rho    ! gradient of the density model
     real(mp),dimension(:)  ,allocatable :: gradient_vp     ! gradient of the velocity model
     real(mp),dimension(:,:),allocatable :: M               ! Mass Matrix
     real(mp),dimension(:,:),allocatable :: Ap
     real(mp),dimension(:,:),allocatable :: App
     real(mp),dimension(:,:),allocatable :: Av

     integer                             :: nb_elem         ! nb of elements
     integer                             :: DoF             ! degree of freedom per elements
     character(len=20)                   :: time_scheme     ! time scheme used
     real(mp)                            :: total_length    ! total length of the domain
     real(mp)                            :: final_time      ! final time of the experiment
     real(mp)                            :: dt              ! time step length
     real(mp)                            :: dx              ! 1D elements size
     real(mp)                            :: alpha           ! penalization
     logical                             :: bernstein       ! logical if T -> Bernstein
     character(len=20)                   :: signal          ! initial perturbation (take flat)
     character(len=20)                   :: boundaries      ! boundary condition
     integer                             :: source_loc      ! source location (at the beginning of an element)
     integer                             :: receiver_loc    ! receiver location (at the beginnig of an element)
     character(len=20)                   :: strategy        ! ATD or DTA
     character(len=20)                   :: scalar_product  ! inner prodct (canonical or M)
     logical                             :: adjoint_test    ! logical for adjoint test
     character(len=20)                   :: animation       ! string for animation
  end type t_fwi

  public  :: t_fwi,init_fwi,free_fwi,one_fwi_step,                      &
             gradient_ATD_computation,gradient_DTA_computation,         &
             update_model,cost_function
  private :: fwi_extract_g

contains

  ! Initializes the fwi type
  subroutine init_fwi(fwi,nb_iter,velocity_ini,density_ini,data_P,data_U,       &
                      nb_elem,DoF,time_scheme,total_length,final_time,alpha,    &
                      bernstein,source_loc,receiver_loc,strategy,               &
                      scalar_product,animation,adjoint_test)
    type(t_fwi)            ,intent(inout) :: fwi
    integer                ,intent(in)    :: nb_iter
    real(mp),dimension(:)  ,intent(in)    :: velocity_ini
    real(mp),dimension(:)  ,intent(in)    :: density_ini
    real(mp),dimension(:,:),intent(in)    :: data_P
    real(mp),dimension(:,:),intent(in)    :: data_U
    integer                ,intent(in)    :: nb_elem
    integer                ,intent(in)    :: DoF
    character(len=*)       ,intent(in)    :: time_scheme
    real(mp)               ,intent(in)    :: total_length
    real(mp)               ,intent(in)    :: final_time
    real(mp)               ,intent(in)    :: alpha
    logical                ,intent(in)    :: bernstein
    integer                ,intent(in)    :: source_loc
    integer                ,intent(in)    :: receiver_loc
    character(len=*)       ,intent(in)    :: strategy
    character(len=*)       ,intent(in)    :: scalar_product
    character(len=*)       ,intent(in)    :: animation
    logical                ,intent(in)    :: adjoint_test

    integer           :: i,j,dv,dp
    real(mp)          :: epsilon_vp
    character(len=50) :: fichier

    fwi%nb_iter=nb_iter

    allocate(fwi%velocity_model(nb_elem))
    allocate(fwi%density_model(nb_elem))
    allocate(fwi%gradient_rho(nb_elem))
    allocate(fwi%gradient_vp(nb_elem))

    dv=max(nb_elem/size(velocity_ini),1)
    do i=1,size(velocity_ini)
       do j=(i-1)*dv+1,i*dv
          fwi%velocity_model(j)=velocity_ini(i)
       end do
    end do
    do j=size(velocity_ini)*dv+1,nb_elem
       fwi%velocity_model(j)=velocity_ini(size(velocity_ini))
    end do

    dp=max(nb_elem/size(density_ini),1)
    do i=1,size(density_ini)
       do j=(i-1)*dp+1,i*dp
          fwi%density_model(j)=density_ini(i)
       end do
    end do
    do j=size(density_ini)*dp+1,nb_elem
       fwi%density_model(j)=density_ini(size(density_ini))
    end do

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
    fwi%source_loc=source_loc
    fwi%receiver_loc=receiver_loc
    fwi%strategy=strategy
    fwi%scalar_product=scalar_product
    fwi%animation=animation
    fwi%adjoint_test=adjoint_test

    fwi%signal='flat'
    fwi%boundaries='ABC'

    write(fichier,"(A,A,I0,'.dat')") "../Files/",'VP',0
    open(unit=28,file=fichier)
    do i=1,fwi%nb_elem
       write(28,*) i,fwi%velocity_model(i)
    end do

  end subroutine init_fwi

  ! Deallocates all fwi type allocatable arrays
  subroutine free_fwi(fwi)
    type(t_fwi),intent(inout) :: fwi

    deallocate(fwi%velocity_model)
    deallocate(fwi%density_model)
    deallocate(fwi%data_P)
    deallocate(fwi%data_U)
    deallocate(fwi%gradient_vp)
    deallocate(fwi%gradient_rho)
  end subroutine free_fwi


  ! Performs one fwi iteration
  subroutine one_fwi_step(fwi)
    type(t_fwi),intent(inout) :: fwi
    integer                   :: current_time_step
    real(mp)                  :: t
    real(mp)                  :: t1,t2

    real(mp),dimension(fwi%DoF*fwi%nb_elem) :: RTM
    integer                                 :: i,j
    type(t_adjoint_test)                    :: test


    call init_acoustic_problem(fwi%forward,fwi%nb_elem,fwi%DoF,fwi%time_scheme, &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%source_loc,      &
         fwi%receiver_loc)

    fwi%n_time_step=fwi%forward%n_time_step
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
    fwi%P_received(:,:)=0.0_mp
    fwi%U_received(:,:)=0.0_mp

    if (fwi%adjoint_test) then
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
       fwi%P_received(current_time_step,2)=                                     &
                                    fwi%forward%P(fwi%DoF*(fwi%receiver_loc-1)+1)

       fwi%U_received(current_time_step,1)=t
       fwi%U_received(current_time_step,2)=                                     &
                                    fwi%forward%U(fwi%DoF*(fwi%receiver_loc-1)+1)

       if (fwi%adjoint_test) then
          fwi%FP(:,current_time_step)=fwi%forward%RHSp
          fwi%FU(:,current_time_step)=fwi%forward%RHSu
          if (fwi%time_scheme.eq.'RK4') then
             fwi%FP_half(:,current_time_step)=fwi%forward%RHSp_half
             fwi%FU_half(:,current_time_step)=fwi%forward%RHSu_half
          end if
       end if
    end do

    if (fwi%adjoint_test) then
       fwi%FP(:,fwi%forward%n_time_step)=0.0_mp
       fwi%FU(:,fwi%forward%n_time_step)=0.0_mp
       fwi%FP(:,0)=0.0_mp
       fwi%FU(:,0)=0.0_mp
       call fwi_extract_g(fwi)
    end if

    if (fwi%animation.eq.'fwi_forward') then
       do i=0,fwi%n_time_step
          if (modulo(i,10).eq.0) then

             call print_vect(fwi%P(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,     &
                             fwi%bernstein,i,'P')

          else if (i.eq.0) then

             call print_vect(fwi%P(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,     &
                             fwi%bernstein,i,'P')

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

    call init_adjoint_problem(fwi%backward,fwi%nb_elem,fwi%DoF,fwi%time_scheme, &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%source_loc,      &
         fwi%receiver_loc,fwi%data_P,fwi%data_U,fwi%P_received,fwi%U_received,  &
         fwi%strategy,fwi%scalar_product)

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
       t=current_time_step*fwi%backward%dt
       call one_backward_step(fwi%backward,t)
       fwi%QP(:,current_time_step)=fwi%backward%P
       fwi%QU(:,current_time_step)=fwi%backward%U


       if (fwi%adjoint_test) then
          if (fwi%time_scheme.eq.'RK4') then
             fwi%DP(:,current_time_step)=fwi%backward%RHSp
             fwi%DU(:,current_time_step)=fwi%backward%RHSu
          else if (fwi%time_scheme.eq.'AB3') then
             fwi%DP(:,current_time_step)=fwi%backward%RHSp
             fwi%DU(:,current_time_step)=fwi%backward%RHSu
          end if
       end if
    end do

    if (fwi%adjoint_test) then
       fwi%DP(:,0)=0.0_mp
       fwi%DU(:,0)=0.0_mp
       fwi%DP(:,fwi%backward%n_time_step)=0.0_mp
       fwi%DU(:,fwi%backward%n_time_step)=0.0_mp
    end if

    if (fwi%animation.eq.'fwi_backward') then
       do i=0,fwi%n_time_step
          if (modulo(i,10).eq.0) then
             call print_vect(fwi%QP(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,    &
                             fwi%bernstein,i,'QP')
          else if (i.eq.0) then
             call print_vect(fwi%QP(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,    &
                             fwi%bernstein,i,'QP')
          end if
       end do
    end  if

    if (fwi%adjoint_test) then
       print*,'%%%%%%% ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
       print*,'Inner Product PU/DPDU',inner_product(fwi%P,fwi%U,fwi%DP,fwi%DU)
       print*,'Inner Product QPQU/GPGU',inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU)
       print*,'Difference relative :',(inner_product(fwi%P,fwi%U,fwi%DP,fwi%DU) &
            -inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU))                        &
            /inner_product(fwi%QP,fwi%QU,fwi%GP,fwi%GU)

       print*,'test',size(fwi%P,1),size(fwi%P,2)

       print*,'M Inner Product PU/DPDU',inner_product_M(fwi%P,fwi%U,   &
                                                        fwi%DP,fwi%DU, &
                                                        test%M)
       print*,'M Inner Product QPQU/GPGU',inner_product_M(fwi%QP,fwi%QU, &
                                                          fwi%GP,fwi%GU, &
                                                          test%M)
       
       print*,'Difference relative :',(inner_product_M(fwi%P,fwi%U,   &
                                                       fwi%DP,fwi%DU, &
                                                       test%M)        &
                                      -inner_product_M(fwi%QP,fwi%QU, &
                                                       fwi%GP,fwi%GU, &
                                                       test%M))       &
                                      /inner_product_M(fwi%QP,fwi%QU, &
                                                       fwi%GP,fwi%GU, &
                                                       test%M)
       print*,'%%%%%% END ADJOINT TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    end if

    call free_adjoint_problem(fwi%backward)

    if (fwi%strategy.eq.'ATD') then
       call gradient_ATD_computation(fwi)
    else if (fwi%strategy.eq.'DTA') then
       call gradient_DTA_computation(fwi)
    end if
    call update_model(fwi)

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
       deallocate(fwi%DP)
       deallocate(fwi%DU)
    end if

  end subroutine one_fwi_step


  ! Evaluates the combined RHS for a specific time scheme (EF=G)
  subroutine fwi_extract_g(fwi)
    type(t_fwi),intent(inout) :: fwi

    real(mp),dimension(size(fwi%P,1)) :: Pk1,Pk2
    real(mp),dimension(size(fwi%P,1)) :: Uk1,Uk2
    integer                           :: i
    real(mp),dimension(size(fwi%P,1)) :: zero
    zero=0.0_mp
    fwi%GP=0.0_mp
    fwi%GU=0.0_mp

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

       Pk1=0.0_mp
       Pk2=0.0_mp
       Uk1=0.0_mp
       Uk2=0.0_mp

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


  ! Evaluates the FWI gradient for the ATD strategy
  subroutine gradient_ATD_computation(fwi)
    type(t_fwi),intent(inout)               :: fwi
    real(mp),dimension(fwi%DoF)             :: dtu,dtp,qu,qp
    real(mp),dimension(fwi%DoF*fwi%nb_elem) :: u0,u1,u2
    real(mp),dimension(fwi%DoF*fwi%nb_elem) :: p0,p1,p2
    real(mp),dimension(fwi%DoF*fwi%nb_elem) :: dtp_glob

    integer  :: nb_elem,DoF,beg_node,end_node
    integer  :: i,j,k
    real(mp) :: moyenne


    nb_elem=fwi%nb_elem   ! For sake of
    DoF=fwi%DoF           ! readibility

    fwi%gradient_vp=0.0_mp
    fwi%gradient_rho=0.0_mp

    do j=0,fwi%n_time_step
       do i=1,nb_elem
          if (j.eq.0) then
             beg_node=(i-1)*DoF+1
             end_node=i*DoF

             dtu=(fwi%U(beg_node:end_node,j+1)-fwi%U(beg_node:end_node,j))      &
                  /(fwi%dt)

             dtp=(fwi%P(beg_node:end_node,j+1)-fwi%P(beg_node:end_node,j))      &
                  /(fwi%dt)

             qp=fwi%QP(beg_node:end_node,j)
             qu=fwi%QU(beg_node:end_node,j)

          else if (j.eq.fwi%n_time_step) then
             beg_node=(i-1)*DoF+1
             end_node=i*DoF

             dtu=(fwi%U(beg_node:end_node,j)-fwi%U(beg_node:end_node,j-1))      &
                  /(fwi%dt)

             dtp=(fwi%P(beg_node:end_node,j)-fwi%P(beg_node:end_node,j-1))      &
                  /(fwi%dt)

             qp=fwi%QP(beg_node:end_node,j)
             qu=fwi%QU(beg_node:end_node,j)
          else
             beg_node=(i-1)*DoF+1
             end_node=i*DoF

             dtu=(fwi%U(beg_node:end_node,j+1)-fwi%U(beg_node:end_node,j-1))    &
                  /(2*fwi%dt)

             dtp=(fwi%P(beg_node:end_node,j+1)-fwi%P(beg_node:end_node,j-1))    &
                  /(2*fwi%dt)

             qp=fwi%QP(beg_node:end_node,j)
             qu=fwi%QU(beg_node:end_node,j)
          end if


          if (fwi%bernstein) then

             dtp=matmul(B2L,dtp)
             dtu=matmul(B2L,dtu)
             qp=matmul(B2L,qp)
             qu=matmul(B2L,qu)

          end if

          moyenne=0.0_mp
          do k=1,DoF
             moyenne=dtp(k)*qp(k)
          end do
          moyenne=moyenne/DoF

          fwi%gradient_vp(i)=fwi%gradient_vp(i)                              &
               -2.0_mp/(fwi%density_model(i)*fwi%velocity_model(i)**3)*moyenne

       end do
    end do

       fwi%gradient_vp=fwi%gradient_vp*fwi%dx*fwi%dt

     end subroutine gradient_ATD_computation


     ! Evaluate the FWI gradient for the DTA strategy
     subroutine gradient_DTA_computation(fwi)
       type(t_fwi),intent(inout)               :: fwi
       real(mp),dimension(fwi%DoF)             :: dtu,dtp,qu,qp
       real(mp),dimension(fwi%DoF*fwi%nb_elem) :: u0,u1,u2
       real(mp),dimension(fwi%DoF*fwi%nb_elem) :: p0,p1,p2
       real(mp),dimension(fwi%DoF*fwi%nb_elem) :: dtp_glob

       integer                                         :: nb_elem
       integer                                         :: DoF
       integer                                         :: beg_node
       integer                                         :: end_node
       integer                                         :: i,j,k
       real(mp)                                        :: one_parameter
       real(mp)                                        :: moyenne
       real(mp),dimension(fwi%DoF,fwi%DoF)             :: M_loc
       real(mp),dimension(fwi%DoF,fwi%DoF*fwi%nb_elem) :: Ap_loc
       real(mp),dimension(fwi%DoF,fwi%DoF*fwi%nb_elem) :: App_loc
       real(mp),dimension(fwi%DoF,fwi%DoF*fwi%nb_elem) :: Av_loc

       type(sparse_matrix)                             :: Ap_sparse_loc
       type(sparse_matrix)                             :: App_sparse_loc
       type(sparse_matrix)                             :: Av_sparse_loc


       nb_elem=fwi%nb_elem   ! For sake of
       DoF=fwi%DoF           ! readibility
       fwi%gradient_vp=0.0_mp
       fwi%gradient_rho=0.0_mp


       if (fwi%time_scheme.eq.'RK4') then
          do j=0,fwi%n_time_step
             do i=1,nb_elem
                M_loc=fwi%M((i-1)*DoF+1:i*DoF,(i-1)*DoF+1:i*DoF)

             if (j.eq.0) then
                beg_node=(i-1)*DoF+1
                end_node=i*DoF

                dtu=(fwi%U(beg_node:end_node,j+1)-fwi%U(beg_node:end_node,j))   &
                     /(fwi%dt)

                dtp=(fwi%P(beg_node:end_node,j+1)-fwi%P(beg_node:end_node,j))   &
                     /(fwi%dt)
                
                qp=fwi%QP(beg_node:end_node,j)
                qu=fwi%QU(beg_node:end_node,j)

             else if (j.eq.fwi%n_time_step) then
                beg_node=(i-1)*DoF+1
                end_node=i*DoF

                dtu=(fwi%U(beg_node:end_node,j)-fwi%U(beg_node:end_node,j-1))   &
                     /(fwi%dt)

                dtp=(fwi%P(beg_node:end_node,j)-fwi%P(beg_node:end_node,j-1))   &
                     /(fwi%dt)

                qp=fwi%QP(beg_node:end_node,j)
                qu=fwi%QU(beg_node:end_node,j)
             else
                beg_node=(i-1)*DoF+1
                end_node=i*DoF

                dtu=(fwi%U(beg_node:end_node,j+1)-fwi%U(beg_node:end_node,j-1)) &
                     /(2*fwi%dt)

                dtp=(fwi%P(beg_node:end_node,j+1)-fwi%P(beg_node:end_node,j-1)) &
                     /(2*fwi%dt)

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

             moyenne=0.0_mp
             do k=1,DoF
                moyenne=dtp(k)*qp(k)
             end do
             moyenne=moyenne/DoF

             fwi%gradient_vp(i)=fwi%gradient_vp(i)                              &
                  -2.0_mp/(fwi%density_model(i)*fwi%velocity_model(i)**3)*moyenne
          end do
       end do

    else

       fwi%gradient_vp=fwi%gradient_vp*fwi%dx*fwi%dt

       if ((fwi%strategy.eq.'DTA').and.(fwi%time_scheme.eq.'AB3')) then
          fwi%gradient_vp=0.0_mp

          fwi%App(1,1)=fwi%App(1,1)/fwi%velocity_model(1)
          fwi%App(fwi%DoF*fwi%nb_elem,fwi%DoF*fwi%nb_elem)=                     &
               fwi%App(fwi%DoF*fwi%nb_elem,fwi%DoF*fwi%nb_elem)/                &
               fwi%velocity_model(fwi%DoF)

             do i=1,nb_elem
                Ap_loc=fwi%Ap((i-1)*DoF+1:i*DoF,:)/(0.5_mp*fwi%velocity_model(i))
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
                   u1=0.0_mp
                   u2=0.0_mp
                   p1=0.0_mp
                   p2=0.0_mp
                else if (j.eq.1) then
                   u1=(fwi%U(:,j-1))
                   p1=(fwi%P(:,j-1))
                   u2=0.0_mp
                   p2=0.0_mp
                else
                   u1=(fwi%U(:,j-1))
                   p1=(fwi%P(:,j-1))
                   u2=(fwi%U(:,j-2))
                   p2=(fwi%P(:,j-2))

                   dtp=-23.0_mp/12.0_mp*(sparse_matmul(App_sparse_loc,p0)   &
                                   +sparse_matmul(Ap_sparse_loc,u0))  &
                        +16.0_mp/12.0_mp*(sparse_matmul(App_sparse_loc,p1)  &
                                    +sparse_matmul(Ap_sparse_loc,u1)) &
                        -5.0_mp/12.0_mp*(sparse_matmul(App_sparse_loc,p2)   &
                                    +sparse_matmul(Ap_sparse_loc,u2))

                   moyenne=0.0_mp
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

  ! Updates the model thanks to the gradient
  subroutine update_model(fwi)
    type(t_fwi),intent(inout) :: fwi

    integer                         :: i,j
    integer,parameter               :: nb_section=9
    integer,dimension(nb_section,2) :: lim_section
    real(mp),dimension(nb_section)  :: grad
    real(mp),parameter              :: beta=0.03_mp
    character(len=50)               :: fichier


    grad=0.0_mp

    do i=1,nb_section
       lim_section(i,1)=(i-1)*int(fwi%nb_elem/nb_section)+1
       lim_section(i,2)=i*int(fwi%nb_elem/nb_section)
    end do
    lim_section(nb_section,2)=fwi%nb_elem

    grad(1:3)=0.0_mp
    do i=4,nb_section
       grad(i)=sum(fwi%gradient_vp(lim_section(i,1):lim_section(i,2)))
    end do

    if (maxval(abs(grad)).eq.0.0_mp) then
       print*,'The gradient is null, end of the fwi algorithm'
       STOP
    else
       grad=(1.0_mp/maxval(abs(grad)))*grad
       !grad=(1.0_mp/grad(50))*grad
    end if

    open(unit=19,file='grad.dat')
    do i=1,size(grad)
       write(19,*) i,grad(i)
    end do

    do i=1,nb_section
       do j=lim_section(i,1),lim_section(i,2)
          fwi%velocity_model(j)=fwi%velocity_model(j)-beta*grad(i)
          if (fwi%velocity_model(j).lt.1.0_mp) then
             fwi%velocity_model(j)=1.0_mp
          else if (fwi%velocity_model(j).gt.2.0_mp) then
             fwi%velocity_model(j)=2.0_mp
          end if
       end do
    end do

    write(fichier,"(A,A,I0,'.dat')") "../Files/",'VP',fwi%current_iter
    open(unit=28,file=fichier)
    do i=1,fwi%nb_elem
       write(28,*) i,fwi%velocity_model(i)
    end do
  end subroutine update_model


  ! Evaluates the Cost function
  subroutine cost_function(P_received,data_P,final_time,iter)
    real(mp),dimension(:,:),intent(in) :: P_received
    real(mp),dimension(:,:),intent(in) :: data_P
    real(mp)               ,intent(in) :: final_time
    integer                ,intent(in) :: iter
    real(mp)                           :: CF
    real(mp)                           :: eval_P_received,eval_data_P
    real(mp)                           :: dt,t
    integer                            :: j,i,n_time_step

    dt=P_received(2,1)-P_received(1,1)
    dt=min(dt,data_P(2,1)-data_P(1,1))

    n_time_step=max(size(P_received,1),size(data_P,1))

    CF=0.0_mp

    do j=2,n_time_step-1
       t=(j-1)*dt
       eval_P_received=0.0_mp
       eval_data_P=0.0_mp

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
    CF=0.5_mp*CF

    write(33,*) iter,CF

    if (CF.eq.0.0_mp) then
       print*,'The gradient is null, end of the FWI algorithm'
       STOP
    end if

  end subroutine cost_function

end module m_fwi
