! This module synthtizes the adjoint test
module m_adjoint_test
  use m_kind
  use m_file_function
  use m_matrix
  use m_time_scheme
  implicit none

  ! The adjoint test type
  type t_adjoint_test
     integer                             :: n_time_step      ! number of time step 
     integer                             :: size_v           ! size of the spatial vectors (DoF*nb_elem)
     character(len=20)                   :: time_scheme      ! tiume scheme used
     real(mp),dimension(:,:),allocatable :: P,P2             ! pressure vectors in space and time
     real(mp),dimension(:,:),allocatable :: U,U2             ! velocity vectors in space and time
     real(mp),dimension(:,:),allocatable :: QP,QP2           ! adjoint pressure vectors in space and time
     real(mp),dimension(:,:),allocatable :: QU,QU2           ! adjoint velocity vectors ine space and time
     real(mp),dimension(:,:),allocatable :: FP,FU            ! RHS in space and time
     real(mp),dimension(:,:),allocatable :: FP_half,FU_half  ! RHS in space and time on half time steps (for RK4)
     real(mp),dimension(:,:),allocatable :: DP,DU            ! adjoint RHS in space and time
     real(mp),dimension(:,:),allocatable :: DP_half,DU_half  ! adjoint RHS in space and time on half time steps (for RK4)
     real(mp),dimension(:,:),allocatable :: GP,GU            ! Combined forward RHS (G=EF)
     real(mp)                            :: dt               ! time step length
     type(sparse_matrix)                 :: Ap,Ap_star
     type(sparse_matrix)                 :: Av,Av_star
     type(sparse_matrix)                 :: App,App_star
     type(sparse_matrix)                 :: M                ! Mass matrix 
     type(sparse_matrix)                 :: Minv             ! Inverse Mass matrix
     integer                             :: nb_elem          ! number of elements
     integer                             :: DoF              ! Degree of freedom
     real(mp)                            :: dx               ! step space length
  end type t_adjoint_test


  contains


    ! Initalizes the adjoint test type
    subroutine init_adjoint_test(test,n_time_step,time_scheme,dt,Ap,Av,App,     &
                                 M,Minv,nb_elem,DoF,dx,FP,FU,FP_half,FU_half,DP,&
                                 DU)

      type(t_adjoint_test),intent(inout) :: test
      integer             ,intent(in)    :: n_time_step
      character(len=*)    ,intent(in)    :: time_scheme
      real(mp)            ,intent(in)    :: dt
      type(sparse_matrix) ,intent(in)    :: Ap
      type(sparse_matrix) ,intent(in)    :: Av
      type(sparse_matrix) ,intent(in)    :: App
      type(sparse_matrix) ,intent(in)    :: M
      type(sparse_matrix) ,intent(in)    :: Minv
      integer             ,intent(in)    :: nb_elem
      integer             ,intent(in)    :: DoF
      real(mp)            ,intent(in)    :: dx

      real(mp),dimension(:,:) ,optional      :: FP,FU
      real(mp),dimension(:,:) ,optional      :: FP_half,FU_half
      real(mp),dimension(:,:) ,optional      :: DP,DU

      integer :: i,j
      type(sparse_matrix) :: Id

      test%size_v=size(Ap%IA)-1
      test%n_time_step=n_time_step
      test%time_scheme=time_scheme
      test%dt=dt
      test%Ap=Ap
      call transpose_sparse(Ap,test%Ap_star)
      test%App=App
      call transpose_sparse(App,test%App_star)
      test%Av=Av
      call transpose_sparse(Av,test%Av_star)
      test%nb_elem=nb_elem
      test%DoF=DoF
      test%dx=dx
      test%M=M
      test%Minv=Minv

      test%App_star=sparse_matmul(test%App_star,test%M)
      test%Ap_star=sparse_matmul(test%Ap_star,test%M)
      test%Av_star=sparse_matmul(test%Av_star,test%M)

      test%App_star=sparse_matmul(test%Minv,test%App_star)
      test%Ap_star=sparse_matmul(test%Minv,test%Ap_star)
      test%Av_star=sparse_matmul(test%Minv,test%Av_star)


      allocate(test%P(0:n_time_step,test%size_v),test%P2(0:n_time_step,test%size_v))
      allocate(test%U(0:n_time_step,test%size_v),test%U2(0:n_time_step,test%size_v))
      allocate(test%QP(0:n_time_step,test%size_v),test%QP2(0:n_time_step,test%size_v))
      allocate(test%QU(0:n_time_step,test%size_v),test%QU2(0:n_time_step,test%size_v))
      allocate(test%FP(0:n_time_step,test%size_v),test%FU(0:n_time_step,test%size_v))
      allocate(test%FP_half(1:n_time_step,test%size_v),test%FU_half(1:n_time_step,test%size_v))
      allocate(test%DP(0:n_time_step,test%size_v),test%DU(0:n_time_step,test%size_v))
      allocate(test%DP_half(1:n_time_step,test%size_v),test%DU_half(1:n_time_step,test%size_v))
      allocate(test%GP(0:n_time_step,test%size_v),test%GU(0:n_time_step,test%size_v))
      ! call random_number(test%P)
      ! call random_number(test%U)
      ! call random_number(test%QP)
      ! call random_number(test%QU)

      if (present(FP)) then

         do i=0,n_time_step
            do j=1,test%size_v
               test%FP(i,j)=FP(j,i+1)
               test%FU(i,j)=FU(j,i+1)
               test%DP(i,j)=DP(j,i+1)
               test%DU(i,j)=DU(j,i+1)
            end do
         end do

         do i=1,n_time_step
            do j=1,test%size_v
               test%FP_half(i,j)=FP_half(j,i)
               test%FU_half(i,j)=FU_half(j,i)
            end do
         end do

         ! test%FP=0.0_mp
         ! test%FU=0.0_mp
         ! test%FP_half=0.0_mp
         ! test%FU_half=0.0_mp

         ! test%DP=0.0_mp
         ! test%DU=0.0_mp
         ! test%DP_half=0.0_mp
         ! test%DU_half=0.0_mp

         ! call random_number(test%FP(400,100))
         ! call random_number(test%FP_half(400,100))
         ! call random_number(test%FU(400,100))
         ! call random_number(test%FU_half(:,100))
         ! call random_number(test%DP(n_time_step-400,100))
         ! call random_number(test%DP_half(n_time_step-400,100))
         ! call random_number(test%DU(n_time_step-400,100))
         ! call random_number(test%DU_half(n_time_step-400,100))
         
         test%P=0.0_mp
         test%U=0.0_mp
         test%QP=0.0_mp
         test%QU=0.0_mp

         if (test%time_scheme.eq.'RK4') then
            test%FP(0,:)=0.0_mp
            test%FU(0,:)=0.0_mp
            test%DP(n_time_step,:)=0.0_mp
            test%DU(n_time_step,:)=0.0_mp
         else if (test%time_scheme.eq.'AB3') then
            test%FP(0:1,:)=0.0_mp
            test%FU(0:1,:)=0.0_mp
            test%DP(n_time_step:n_time_step-1,:)=0.0_mp
            test%DU(n_time_step:n_time_step-1,:)=0.0_mp

            test%FP(0,:)=0.0_mp
            test%FU(0,:)=0.0_mp
            test%DP(n_time_step,:)=0.0_mp
            test%DU(n_time_step,:)=0.0_mp
         end if


      else
         test%P=0.0_mp
         test%U=0.0_mp
         test%QP=0.0_mp
         test%QU=0.0_mp
         call random_number(test%FP)
         call random_number(test%FP_half)
         call random_number(test%FU)
         call random_number(test%FU_half)
         call random_number(test%DP)
         call random_number(test%DP_half)
         call random_number(test%DU)
         call random_number(test%DU_half)

         ! test%FP=0.0_mp
         ! test%FP_half=0.0_mp

         ! open(unit=7,file='fort.22')
         ! read(7,*) test%FU(0,:)
         ! close(7)

         ! print*,'TEST',maxval(test%FU(0,:))

         ! do i=1,test%n_time_step
         !    test%FU(i,:)=sin(60*i*test%dt)*test%FU(0,:)
         !    test%FU_half(i-1,:)=sin(60*(i-0.5_mp)*test%dt)*test%FU(0,:)
         ! end do

         do i=1,test%n_time_step
            do j=1,test%size_v
               test%FU(i,j)=test%FU(i,j)-0.5_mp
               test%FP(i,j)=test%FP(i,j)-0.5_mp
               test%FU_half(i,j)=test%FU_half(i,j)-0.5_mp
               test%FP_half(i,j)=test%FP_half(i,j)-0.5_mp
               test%DU(i,j)=test%DU(i,j)-0.5_mp
               test%DP(i,j)=test%DP(i,j)-0.5_mp
               test%DU_half(i,j)=test%DU_half(i,j)-0.5_mp
               test%DP_half(i,j)=test%DP_half(i,j)-0.5_mp
            end do
         end do

         if (test%time_scheme.eq.'RK4') then
            test%FP(0,:)=0.0_mp
            test%FU(0,:)=0.0_mp
            test%DP(n_time_step,:)=0.0_mp
            test%DU(n_time_step,:)=0.0_mp
         else if (test%time_scheme.eq.'AB3') then
            test%FP(0:1,:)=0.0_mp
            test%FU(0:1,:)=0.0_mp
            test%DP(n_time_step:n_time_step-1,:)=0.0_mp
            test%DU(n_time_step:n_time_step-1,:)=0.0_mp

            test%FP(0,:)=0.0_mp
            test%FU(0,:)=0.0_mp
            test%DP(n_time_step,:)=0.0_mp
            test%DU(n_time_step,:)=0.0_mp
         end if


         ! do i=0,n_time_step
         !    test%FP(i,:)=i!real(i/n_time_step)-0.5_mp
         !    test%FU(i,:)=i!real(i/n_time_step)-0.5_mp
         ! end do

         ! do i=1,n_time_step
         !    test%FP_half(i,:)=i-0.5_mp!real((i-0.5_mp)/n_time_step)-0.5_mp
         !    test%FU_half(i,:)=i-0.5_mp!real((i-0.5_mp)/n_time_step)-0.5_mp
         ! end do


         ! do i=0,n_time_step
         !    test%DP(i,:)=n_time_step-i!real(i/n_time_step)-0.5_mp
         !    test%DU(i,:)=n_time_step-i!real(i/n_time_step)-0.5_mp
         ! end do

         ! do i=1,n_time_step
         !    test%DP_half(i,:)=n_time_step-i+0.5_mp!real((i-0.5_mp)/n_time_step)-0.5_mp
         !    test%DU_half(i,:)=n_time_step-i+0.5_mp!real((i-0.5_mp)/n_time_step)-0.5_mp
         ! end do
      end if
    end subroutine init_adjoint_test


    ! Proceeds to a forward simulation to initialize P and U in time and space
    subroutine forward_test(test)
      type(t_adjoint_test),intent(inout) :: test


      if (test%time_scheme.eq.'LF') then
         call test_LF_forward(test%P,test%U,test%FP,test%FU_half,test%Ap,       &
                              test%Av,test%App,test%n_time_step,test%size_v,    &
                              test%nb_elem,test%DoF,test%dx)


      else if (test%time_scheme.eq.'RK4') then
         call test_RK4_forward(test%P,test%U,test%FP,test%FU,test%FP_half,      &
                               test%FU_half,test%Ap,test%Av,test%App,           &
                               test%n_time_step,test%size_v,test%nb_elem,       &
                               test%DoF,test%dx)


      else if (test%time_scheme.eq.'AB3') then
         call  test_AB3_forward(test%P,test%U,test%FP,test%FU,test%Ap,test%Av,  &
                                test%App,test%n_time_step,test%size_v,          &
                                test%nb_elem,test%DoF,test%dx)
      else
         print*,'not recongnized time scheme'
      end if
    end subroutine forward_test


    ! Proceeds to a forward simulation to initialize P and U in time and space
    ! using the combined RHS G (G=EF)
    subroutine forward_test2(test)
      type(t_adjoint_test),intent(inout) :: test
      integer :: i

    call extract_g(test)

    if (test%time_scheme.eq.'LF') then
       call test_LF_forward(test%P,test%U,test%FP,test%FU_half,test%Ap,test%Av, &
            test%App,test%n_time_step,test%size_v,test%nb_elem,                 &
            test%DoF,test%dx)


    else if (test%time_scheme.eq.'RK4') then
       call test_RK4_forward2(test%P,test%U,test%GP,test%GU,test%FP_half,       &
            test%FU_half,test%Ap,test%Av,test%App,                              &
            test%n_time_step,test%size_v,test%nb_elem,test%DoF,                 &
            test%dx,test%GP,test%GU)


    else if (test%time_scheme.eq.'AB3') then
       call  test_AB3_forward2(test%P,test%U,test%FP,test%FU,test%Ap,test%Av,   &
            test%App,test%n_time_step,test%size_v,                              &
            test%nb_elem,test%DoF,test%dx,test%GP,test%GU)
    else
       print*,'not recongnized time scheme'
    end if
  end subroutine forward_test2


  ! Proceeds to a Backward simulation to initialize the adjoint space time      &
  ! vectors QP and QU
  subroutine backward_test(test)
    type(t_adjoint_test),intent(inout) :: test

    if (test%time_scheme.eq.'LF') then
       call test_LF_backward(test%QP,test%QU,test%DP,test%DU_half,test%Ap,      &
                             test%Av,test%App,test%n_time_step,test%size_v,     &
                             test%nb_elem,test%DoF,test%dx)

    else if (test%time_scheme.eq.'RK4') then
       call test_RK4_backward(test%QP,test%QU,test%DP,test%DU,test%Av_star,     &
                              test%Ap_star,test%App_star,test%n_time_step,      &
                              test%size_v,test%nb_elem,test%DoF,test%dx)

    else if (test%time_scheme.eq.'AB3') then
       call test_AB3_backward(test%QP,test%QU,test%DP,test%DU,test%Av_star,     &
                              test%Ap_star,test%App_star,test%n_time_step,      &
                              test%size_v,test%nb_elem,test%DoF,test%dx,test%GP,&
                              test%GU)

    else
       print*,'not recongnized time scheme'
    end if
  end subroutine backward_test


  ! G=EF
  subroutine extract_g(test)
    type(t_adjoint_test),intent(inout) :: test

    real(mp),dimension(size(test%P,2)) :: Pk1,Pk2
    real(mp),dimension(size(test%P,2)) :: Uk1,Uk2
    real(mp),dimension(size(test%P,2)) :: zero
    integer                            :: i

    zero=0.0_mp
    test%GP=0.0_mp
    test%GU=0.0_mp

    if (test%time_scheme.eq.'LF') then
       print*,'not possible yet for Leap Frog'
       STOP

    else if (test%time_scheme.eq.'RK4') then
       do i=1,test%n_time_step
          call RK4_forward(test%GP(i,:),test%GU(i,:),                           &
               test%Ap,test%Av,test%App,                                        &
               test%FP(i-1,:),test%FP_half(i,:),test%FP(i,:),                   &
               test%FU(i-1,:),test%FU_half(i,:),test%FU(i,:))
       end do

    else if (test%time_scheme.eq.'AB3') then

       Pk1=0.0_mp
       Pk2=0.0_mp
       Uk1=0.0_mp
       Uk2=0.0_mp

       call  AB3_forward(test%GP(0,:),test%GU(0,:),test%Ap,test%Av,test%App,    &
                         test%FP(0,:),test%FU(0,:),                             &
                         zero,zero,                                             &
                         zero,zero,                                             &
                         Pk1,Pk2,Uk1,Uk2)


       call  AB3_forward(test%GP(1,:),test%GU(1,:),test%Ap,test%Av,test%App,    &
                         test%FP(1,:),test%FU(1,:),                             &
                         test%FP(0,:),test%FU(0,:),                             &
                         zero,zero,                                             &
                         Pk1,Pk2,Uk1,Uk2)


       do i=2,test%n_time_step
          call  AB3_forward(test%GP(i,:),test%GU(i,:),test%Ap,test%Av,test%App, &
                            test%FP(i,:),test%FU(i,:),                          &
                            test%FP(i-1,:),test%FU(i-1,:),                      &
                            test%FP(i-2,:),test%FU(i-2,:),                      &
                            Pk1,Pk2,Uk1,Uk2)
       end do

    else
       print*,'not recongnized time scheme'
    end if
  end subroutine extract_g



  ! Updates P and U by using a Leap Frog time scheme
  subroutine test_LF_forward(P,U,FP,FU_half,Ap,Av,App,n_time_step,              &
                             size_v,nb_elem,DoF,dx)

    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: P
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: U
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FP
    real(mp),dimension(1:n_time_step,size_v),intent(inout) :: FU_half
    type(sparse_matrix)                     ,intent(in)    :: Ap,Av,App
    integer                                 ,intent(in)    :: n_time_step
    integer                                 ,intent(in)    :: size_v
    integer                                 ,intent(in)    :: nb_elem
    integer                                 ,intent(in)    :: DoF
    real(mp)                                ,intent(in)    :: dx

    real(mp),dimension(size(U,2)) :: P_current,U_current
    integer                       :: i

    P(0,:)=0.0_mp
    U(0,:)=0.0_mp
    FP(0,:)=0.0_mp

    do i=1,n_time_step
       ! print*,'i',i

       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call LF_forward(P_current,U_current,Ap,Av,App,FP(i-1,:),FU_half(i,:))
       P(i,:)=P_current
       U(i,:)=U_current
    end do
  end subroutine test_LF_forward


  ! Updates QP and QU by using a Leap Frog time scheme
  subroutine test_LF_backward(QP,QU,DP,DU_half,Av_star,Ap_star,App_star,        &
                              n_time_step,size_v,nb_elem,DoF,dx)
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: QP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: QU
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: DP
    real(mp),dimension(1:n_time_step,size_v),intent(inout) :: DU_half
    type(sparse_matrix)                     ,intent(in)    :: Ap_star,Av_star,App_star
    integer                                 ,intent(in)    :: n_time_step
    integer                                 ,intent(in)    :: size_v
    integer                                 ,intent(in)    :: nb_elem
    integer                                 ,intent(in)    :: DoF
    real(mp)                                ,intent(in)    :: dx

    real(mp),dimension(size(QU,2)) :: QP_current,QU_current
    integer                        :: i

    QP(n_time_step,:)=0.0_mp
    QU(n_time_step,:)=0.0_mp
    DP(n_time_step,:)=0.0_mp

    do i=n_time_step-1,0-1
       QP_current=QP(i+1,:)
       QU_current=QU(i+1,:)
       call LF_forward(QP_current,QU_current,Av_star,Ap_star,App_star,DP(i+1,:),&
                       DU_half(i+1,:))
       QP(i,:)=QP_current
       QU(i,:)=QU_current
    end do
  end subroutine test_LF_backward

  ! Updates P and U by using a RK4 time scheme  
  subroutine test_RK4_forward(P,U,FP,FU,FP_half,FU_half,Ap,Av,App,n_time_step,  &
                              size_v,nb_elem,DoF,dx)
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: U
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: P
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FU
    real(mp),dimension(1:n_time_step,size_v),intent(inout) :: FP_half
    real(mp),dimension(1:n_time_step,size_v),intent(inout) :: FU_half
    type(sparse_matrix)                     ,intent(in)    :: Ap,Av,App
    integer                                 ,intent(in)    :: n_time_step
    integer                                 ,intent(in)    :: size_v
    integer                                 ,intent(in)    :: nb_elem
    integer                                 ,intent(in)    :: DoF
    real(mp)                                ,intent(in)    :: dx

    real(mp),dimension(size(U,2)) :: P_current,U_current
    integer                       :: i

    P(0,:)=0.0_mp
    U(0,:)=0.0_mp
    FP(0,:)=0.0_mp
    FU(0,:)=0.0_mp

    do i=1,n_time_step
       ! print*,'i',i

       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call RK4_forward(P_current,U_current,Ap,Av,App, &
            FP(i-1,:),FP_half(i,:),                    &
            FP(i,:),FU(i-1,:),                         &
            FU_half(i,:),FU(i,:),                      &
            0.0_mp*FP(i,:),0.0_mp*FU(i,:))

       P(i,:)=P_current
       U(i,:)=U_current
    end do
  end subroutine test_RK4_forward


  ! Updates P and U by using a RK4 time scheme by using the combined RHS G (G=EF)
  subroutine test_RK4_forward2(P,U,FP,FU,FP_half,FU_half,Ap,Av,App,n_time_step,&
                                size_v,nb_elem,DoF,dx,GP,GU)
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: U
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: P
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FU
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: GP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: GU
    real(mp),dimension(1:n_time_step,size_v),intent(inout) :: FP_half
    real(mp),dimension(1:n_time_step,size_v),intent(inout) :: FU_half
    type(sparse_matrix)                     ,intent(in)    :: Ap,Av,App
    integer                                 ,intent(in)    :: n_time_step
    integer                                 ,intent(in)    :: size_v
    integer                                 ,intent(in)    :: nb_elem
    integer                                 ,intent(in)    :: DoF
    real(mp)                                ,intent(in)    :: dx

    real(mp),dimension(size(U,2)) :: P_current,U_current
    integer                       :: i

    P(0,:)=0.0_mp
    U(0,:)=0.0_mp
    FP(0,:)=0.0_mp
    FU(0,:)=0.0_mp

    do i=1,n_time_step

       P_current=P(i-1,:)
       U_current=U(i-1,:)

       call RK4_forward(P_current,U_current,Ap,Av,App,                          &
                        0.0_mp*FP(i-1,:),0.0_mp*FP_half(i,:),                   &
                        0.0_mp*FP(i,:),0.0_mp*FU(i-1,:),                        &
                        0.0_mp*FU_half(i,:),0.0_mp*FU(i,:),                     &
                        GP(i,:),GU(i,:))     
       P(i,:)=P_current
       U(i,:)=U_current
    end do
  end subroutine test_RK4_forward2



  ! Updates QP and QU by using a RK4 time scheme  
  subroutine test_RK4_backward(QP,QU,DP,DU,Av_star,Ap_star,App_star, &
                               n_time_step,size_v,nb_elem,DoF,dx)
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: QU
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: QP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: DP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: DU
    type(sparse_matrix)                     ,intent(in)    :: Ap_star
    type(sparse_matrix)                     ,intent(in)    :: Av_star
    type(sparse_matrix)                     ,intent(in)    :: App_star
    integer                                 ,intent(in)    :: n_time_step
    integer                                 ,intent(in)    :: size_v
    integer                                 ,intent(in)    :: nb_elem
    integer                                 ,intent(in)    :: DoF
    real(mp)                                ,intent(in)    :: dx

    real(mp),dimension(size_v)     :: zero
    real(mp),dimension(size(QP,2)) :: QP_current,QU_current
    integer                        :: i

    zero=0.0_mp

    QP(n_time_step,:)=0.0_mp
    QU(n_time_step,:)=0.0_mp
    DP(n_time_step,:)=0.0_mp
    DU(n_time_step,:)=0.0_mp

    do i=n_time_step-1,0,-1
       QP_current=QP(i+1,:)
       QU_current=QU(i+1,:)
       call RK4_forward(QP_current,QU_current,                                  &
                        Av_star,Ap_star,App_star,                               &
                        zero,zero,zero,zero,zero,zero,                          &
                        DP(i,:),DU(i,:))
       QP(i,:)=QP_current
       QU(i,:)=QU_current
    end do
  end subroutine test_RK4_backward


  ! Updates P and U by using a AB3 time scheme
  subroutine test_AB3_forward(P,U,FP,FU,Ap,Av,App,n_time_step,size_v,nb_elem,   &
                              DoF,dx)
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: U
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: P
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FU
    
    type(sparse_matrix),intent(in) :: Ap,Av,App
    integer            ,intent(in) :: n_time_step
    integer            ,intent(in) :: size_v
    integer            ,intent(in) :: nb_elem
    integer            ,intent(in) :: DoF
    real(mp)           ,intent(in) :: dx

    real(mp),dimension(size(U,2)) :: P_current,U_current
    real(mp),dimension(size(U,2)) :: Pk1,Pk2
    real(mp),dimension(size(U,2)) :: Uk1,Uk2
    real(mp),dimension(size(U,2)) :: zero
    integer                       :: i


    Pk1=0.0_mp
    Pk2=0.0_mp
    Uk1=0.0_mp
    Uk2=0.0_mp
    zero=0.0_mp

    P_current=P(0,:)
    U_current=U(0,:)
    call  AB3_forward(P_current,U_current,Ap,Av,App,                            &
                      FP(0,:),FU(0,:),                                          &
                      zero,zero,                                                &
                      zero,zero,                                                &
                      Pk1,Pk2,Uk1,Uk2)
    P(1,:)=P_current
    U(1,:)=U_current

    P_current=P(1,:)
    U_current=U(1,:)
    call  AB3_forward(P_current,U_current,Ap,Av,App,                            &
                      FP(1,:),FU(1,:),                                          &
                      FP(0,:),FU(0,:),                                          &
                      zero,zero,                                                &
                      Pk1,Pk2,Uk1,Uk2)
    P(2,:)=P_current
    U(2,:)=U_current


    do i=3,n_time_step
       ! print*,'i',i

       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call  AB3_forward(P_current,U_current,Ap,Av,App,                         &
                            FP(i-1,:),FU(i-1,:),                                &
                            FP(i-2,:),FU(i-2,:),                                &
                            FP(i-3,:),FU(i-3,:),                                &
                            Pk1,Pk2,Uk1,Uk2)
       P(i,:)=P_current
       U(i,:)=U_current
    end do
  end subroutine test_AB3_forward



  ! Updates P and U by using a AB3 time scheme by using the combined RHS G (G=EF)
  subroutine test_AB3_forward2(P,U,FP,FU,Ap,Av,App,n_time_step,size_v,nb_elem,  &
                               DoF,dx,GP,GU)
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: U
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: P
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: FU
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: GP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: GU

    type(sparse_matrix),intent(in)    :: Ap,Av,App
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real(mp)           ,intent(in)    :: dx

    real(mp),dimension(size(U,2)) :: P_current,U_current
    real(mp),dimension(size(U,2)) :: Pk1,Pk2
    real(mp),dimension(size(U,2)) :: Uk1,Uk2
    real(mp),dimension(size(U,2)) :: zero
    integer                       :: i

    Pk1=0.0_mp
    Pk2=0.0_mp
    Uk1=0.0_mp
    Uk2=0.0_mp
    zero=0.0_mp

    P_current=P(0,:)
    U_current=U(0,:)
    call  AB3_forward(P_current,U_current,Ap,Av,App,                            &
                      zero,zero,                                                &
                      zero,zero,                                                &
                      zero,zero,                                                &
                      Pk1,Pk2,Uk1,Uk2,GP(0,:),GU(0,:))
    P(1,:)=P_current
    U(1,:)=U_current

    P_current=P(1,:)
    U_current=U(1,:)
    call  AB3_forward(P_current,U_current,Ap,Av,App,                            &
                      zero,zero,                                                &   
                      zero,zero,                                                &  
                      zero,zero,                                                &
                      Pk1,Pk2,Uk1,Uk2,GP(1,:),GU(1,:))
    P(2,:)=P_current
    U(2,:)=U_current


    do i=2,n_time_step
       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call  AB3_forward(P_current,U_current,Ap,Av,App,                         &
                         zero,zero,                                             &
                         zero,zero,                                             &
                         zero,zero,                                             &
                         Pk1,Pk2,Uk1,Uk2,GP(i-1,:),GU(i-1,:))
       P(i,:)=P_current
       U(i,:)=U_current
    end do
  end subroutine test_AB3_forward2


  ! Updates QP and QU by using a AB3 time scheme  
  subroutine test_AB3_backward(QP,QU,DP,DU,Av_star,Ap_star,App_star,n_time_step,&
                               size_v,nb_elem,DoF,dx,GP,GU)

    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: QU
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: QP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: DP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: DU
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: GP
    real(mp),dimension(0:n_time_step,size_v),intent(inout) :: GU

    type(sparse_matrix),intent(in) :: Ap_star,Av_star,App_star
    integer            ,intent(in) :: n_time_step
    integer            ,intent(in) :: size_v
    integer            ,intent(in) :: nb_elem
    integer            ,intent(in) :: DoF
    real(mp)           ,intent(in) :: dx

    real(mp),dimension(size(QU,2)) :: QP_current,QU_current
    real(mp),dimension(size(QU,2)) :: QPk1,QPk2
    real(mp),dimension(size(QU,2)) :: QUk1,QUk2
    real(mp),dimension(size(QU,2)) :: zero
    integer                        :: i

    QPk1=0.0_mp
    QPk2=0.0_mp
    QUk1=0.0_mp
    QUk2=0.0_mp
    zero=0.0_mp

    QP_current=QP(n_time_step,:)
    QU_current=QU(n_time_step,:)
    call  AB3_forward(QP_current,QU_current,Av_star,Ap_star,App_star,           &
                      zero,zero,                                                &
                      zero,zero,                                                &
                      zero,zero,                                                &
                      QPk1,QPk2,QUk1,QUk2,                                      &
                      DP(n_time_step,:),DU(n_time_step,:))
    QP(n_time_step-1,:)=QP_current
    QU(n_time_step-1,:)=QU_current

    QP_current=QP(n_time_step-1,:)
    QU_current=QU(n_time_step-1,:)
    call  AB3_forward(QP_current,QU_current,Av_star,Ap_star,App_star,           &
                      zero,zero,                                                &
                      zero,zero,                                                &
                      zero,zero,                                                &
                      QPk1,QPk2,QUk1,QUk2,                                      &
                      DP(n_time_step-1,:),DU(n_time_step-1,:))
    QP(n_time_step-2,:)=QP_current
    QU(n_time_step-2,:)=QU_current


    do i=n_time_step-3,0,-1
       QP_current=QP(i+1,:)
       QU_current=QU(i+1,:)
       call  AB3_forward(QP_current,QU_current,Av_star,Ap_star,App_star,        &
                            zero,zero,                                          &
                            zero,zero,                                          &
                            zero,zero,                                          &
                            QPk1,QPk2,QUk1,QUk2,                                &
                            DP(i+1,:),DU(i+1,:) )
       QP(i,:)=QP_current
       QU(i,:)=QU_current
    end do
  end subroutine test_AB3_backward


  ! Evaluates the inner product <P,U / DP,DU>
  function inner_product(P,U,DP,DU)
    real(mp),dimension(:,:),intent(in) :: U,P
    real(mp),dimension(:,:),intent(in) :: DU,DP
    real(mp)                           :: inner_product
    integer                            :: i,j
    integer                            :: n_time_step

    n_time_step=size(U,1)
    inner_product=0.0_mp

    do i=1,n_time_step
       do j=1,size(U,2)
          inner_product=inner_product+P(i,j)*DP(i,j)+U(i,j)*DU(i,j)
       end do
    end do
  end function inner_product


  ! Evaluates the inner product <M (P,U) / DP, DU>
  function inner_product_M(P,U,DP,DU,M)
    real(mp),dimension(:,:),intent(in) :: U,P
    real(mp),dimension(:,:),intent(in) :: DU,DP
    type(sparse_matrix)    ,intent(in) :: M
    real(mp)                           :: inner_product_M
    integer                            :: i,j
    integer                            :: n_time_step
    real(mp),dimension(size(P,1))      :: MP,MU

    n_time_step=size(P,2)
    inner_product_M=0.0

    do i=1,n_time_step

       MP=sparse_matmul(M,P(:,i))
       MU=sparse_matmul(M,U(:,i))

       do j=1,size(P,1)
          inner_product_M=inner_product_M+MP(j)*DP(j,i)+MU(j)*DU(j,i)
       end do
    end do
  end function inner_product_M

end module m_adjoint_test
