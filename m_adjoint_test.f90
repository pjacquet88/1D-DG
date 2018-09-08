module m_adjoint_test
  use m_file_function
  use m_matrix
  use m_time_scheme
  implicit none

  type t_adjoint_test
     integer                         :: n_time_step
     integer                         :: size_v
     character(len=20)               :: time_scheme
     real,dimension(:,:),allocatable :: P,P2
     real,dimension(:,:),allocatable :: U,U2
     real,dimension(:,:),allocatable :: QP,QP2
     real,dimension(:,:),allocatable :: QU,QU2
     real,dimension(:,:),allocatable :: FP,FU
     real,dimension(:,:),allocatable :: FP_half,FU_half
     real,dimension(:,:),allocatable :: DP,DU
     real,dimension(:,:),allocatable :: DP_half,DU_half
     real,dimension(:,:),allocatable :: GP,GU
     real                            :: dt
     type(sparse_matrix)             :: Ap,tAp
     type(sparse_matrix)             :: Av,tAv
     type(sparse_matrix)             :: App,tApp

     integer                         :: nb_elem
     integer                         :: DoF
     real                            :: dx
  end type t_adjoint_test


  contains
  

    subroutine init_adjoint_test(test,n_time_step,time_scheme,dt,Ap,Av,App,       &
         nb_elem,DoF,dx)
    type(t_adjoint_test),intent(inout) :: test
    integer             ,intent(in)    :: n_time_step
    character(len=*)    ,intent(in)    :: time_scheme
    real                ,intent(in)    :: dt
    type(sparse_matrix) ,intent(in)    :: Ap
    type(sparse_matrix) ,intent(in)    :: Av
    type(sparse_matrix) ,intent(in)    :: App
    integer             ,intent(in)    :: nb_elem
    integer             ,intent(in)    :: DoF
    real                ,intent(in)    :: dx

    integer :: i,j
    real,dimension(:,:),allocatable :: Ap_full,tAp_full,ttAp_full

    test%size_v=size(Ap%IA)-1
    test%n_time_step=n_time_step
    test%time_scheme=time_scheme
    test%dt=dt
    test%Ap=Ap
    call transpose_sparse(Ap,test%tAp)
    test%App=App
    call transpose_sparse(App,test%tApp)
    test%Av=Av
    call transpose_sparse(Av,test%tAv)
    test%nb_elem=nb_elem
    test%DoF=DoF
    test%dx=dx


    
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
    test%P=0.0
    test%U=0.0
    test%QP=0.0
    test%QU=0.0
    call random_number(test%FP)
    call random_number(test%FP_half)
    call random_number(test%FU)
    call random_number(test%FU_half)
    call random_number(test%DP)
    call random_number(test%DP_half)
    call random_number(test%DU)
    call random_number(test%DU_half)

    ! test%FP=0.0
    ! test%FP_half=0.0

    ! open(unit=7,file='fort.22')
    ! read(7,*) test%FU(0,:)
    ! close(7)

    ! print*,'TEST',maxval(test%FU(0,:))

    ! do i=1,test%n_time_step
    !    test%FU(i,:)=sin(60*i*test%dt)*test%FU(0,:)
    !    test%FU_half(i-1,:)=sin(60*(i-0.5)*test%dt)*test%FU(0,:)
    ! end do
    
    do i=1,test%n_time_step
       do j=1,test%size_v
          test%FU(i,j)=test%FU(i,j)-0.5
          test%FP(i,j)=test%FP(i,j)-0.5
          test%FU_half(i,j)=test%FU_half(i,j)-0.5
          test%FP_half(i,j)=test%FP_half(i,j)-0.5
          test%DU(i,j)=test%DU(i,j)-0.5
          test%DP(i,j)=test%DP(i,j)-0.5
          test%DU_half(i,j)=test%DU_half(i,j)-0.5
          test%DP_half(i,j)=test%DP_half(i,j)-0.5
       end do
    end do

    do i=0,n_time_step
       test%FP(i,:)=i!real(i/n_time_step)-0.5
       test%FU(i,:)=i!real(i/n_time_step)-0.5
    end do

    do i=1,n_time_step
       test%FP_half(i,:)=i-0.5!real((i-0.5)/n_time_step)-0.5
       test%FU_half(i,:)=i-0.5!real((i-0.5)/n_time_step)-0.5
    end do

    
    do i=0,n_time_step
       test%DP(i,:)=n_time_step-i!real(i/n_time_step)-0.5
       test%DU(i,:)=n_time_step-i!real(i/n_time_step)-0.5
    end do

    do i=1,n_time_step
       test%DP_half(i,:)=n_time_step-i+0.5!real((i-0.5)/n_time_step)-0.5
       test%DU_half(i,:)=n_time_step-i+0.5!real((i-0.5)/n_time_step)-0.5
    end do

  end subroutine init_adjoint_test

  
  subroutine forward_test(test)
    type(t_adjoint_test),intent(inout) :: test


    if (test%time_scheme.eq.'LF') then
       call test_LF_forward(test%P,test%U,test%FP,test%FU_half,test%Ap,test%Av, &
            test%App,test%n_time_step,test%size_v,test%nb_elem, &
            test%DoF,test%dx)


    else if (test%time_scheme.eq.'RK4') then
       call test_RK4_forward(test%P,test%U,test%FP,test%FU,test%FP_half,        &
            test%FU_half,test%Ap,test%Av,test%App,             &
            test%n_time_step,test%size_v,test%nb_elem,test%DoF,&
            test%dx)


    else if (test%time_scheme.eq.'AB3') then
       call  test_AB3_forward(test%P,test%U,test%FP,test%FU,test%Ap,test%Av,    &
            test%App,test%n_time_step,test%size_v,            &
            test%nb_elem,test%DoF,test%dx)

    else
       print*,'not recongnized time scheme'
    end if

  end subroutine forward_test

    subroutine forward_test2(test)
    type(t_adjoint_test),intent(inout) :: test

    
    call extract_g(test)

    if (test%time_scheme.eq.'LF') then
       call test_LF_forward(test%P,test%U,test%FP,test%FU_half,test%Ap,test%Av, &
            test%App,test%n_time_step,test%size_v,test%nb_elem, &
            test%DoF,test%dx)


    else if (test%time_scheme.eq.'RK4') then
       call test_RK4_forward2(test%P,test%U,test%GP,test%GU,test%FP_half,        &
            test%FU_half,test%Ap,test%Av,test%App,             &
            test%n_time_step,test%size_v,test%nb_elem,test%DoF,&
            test%dx,test%GP,test%GU)


    else if (test%time_scheme.eq.'AB3') then
       call  test_AB3_forward(test%P,test%U,test%FP,test%FU,test%Ap,test%Av,    &
            test%App,test%n_time_step,test%size_v,            &
            test%nb_elem,test%DoF,test%dx)

    else
       print*,'not recongnized time scheme'
    end if

  end subroutine forward_test2
     
  subroutine backward_test(test)
    type(t_adjoint_test),intent(inout) :: test

    if (test%time_scheme.eq.'LF') then
       call test_LF_backward(test%QP,test%QU,test%DP,test%DU_half,test%Ap,test%Av, &
            test%App,test%n_time_step,test%size_v,test%nb_elem, &
            test%DoF,test%dx)
    
    else if (test%time_scheme.eq.'RK4') then
       call test_RK4_backward(test%QP,test%QU,test%DP,test%DU,test%tAv,test%tAp,  &
            test%tApp,test%n_time_step,test%size_v,test%nb_elem,test%DoF,test%dx)

    else if (test%time_scheme.eq.'AB3') then
       call test_AB3_backward(test%QP,test%QU,test%DP,test%DU,test%tAv,test%tAp,&
            test%tApp,test%n_time_step,test%size_v,            &
            test%nb_elem,test%DoF,test%dx)

    else
       print*,'not recongnized time scheme'
    end if

  end subroutine backward_test

  subroutine extract_g(test)
    type(t_adjoint_test),intent(inout) :: test

    integer :: i

    test%GP=0.0
    test%GU=0.0
    
    if (test%time_scheme.eq.'LF') then
       print*,'not possible yet for Leap Frog'
       STOP

       
    else if (test%time_scheme.eq.'RK4') then
       do i=1,test%n_time_step
          call RK4_forward(test%GP(i,:),test%GU(i,:),                        &
               test%Ap,test%Av,test%App,                                     &
               test%FP(i-1,:),test%FP_half(i,:),test%FP(i,:),                &
               test%FU(i-1,:),test%FU_half(i,:),test%FU(i,:))
       end do

    else if (test%time_scheme.eq.'AB3') then
       do i=1,test%n_time_step
          call  AB3_forward(test%GP(i,:),test%GU(i,:),test%Ap,test%Av,test%App, &
               test%FP(i,:),test%FU(i,:),                                    &
               test%FP(i-1,:),test%FP(i-2,:),test%FU(i-1,:),test%FU(i-2,:))
       end do
    else
       print*,'not recongnized time scheme'
    end if

  end subroutine extract_g

  subroutine test_LF_forward(P,U,FP,FU_half,Ap,Av,App,n_time_step,       &
       size_v,nb_elem,DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: P
    real,dimension(0:n_time_step,size_v),intent(inout) :: U
    real,dimension(0:n_time_step,size_v),intent(inout) :: FP
    real,dimension(1:n_time_step,size_v),intent(inout) :: FU_half
    type(sparse_matrix),intent(in)    :: Ap,Av,App
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(U,2)) :: P_current,U_current
    integer                   :: i

    P(0,:)=0.0
    U(0,:)=0.0
    FP(0,:)=0.0
    
    do i=1,n_time_step
       ! print*,'i',i

       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call LF_forward(P_current,U_current,Ap,Av,App,FP(i-1,:),FU_half(i,:))
       P(i,:)=P_current
       U(i,:)=U_current
    end do
  end subroutine test_LF_forward

  subroutine test_LF_backward(QP,QU,DP,DU_half,tAv,tAp,tApp,n_time_step,      &
       size_v,nb_elem,DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: QP
    real,dimension(0:n_time_step,size_v),intent(inout) :: QU
    real,dimension(0:n_time_step,size_v),intent(inout) :: DP
    real,dimension(1:n_time_step,size_v),intent(inout) :: DU_half
    type(sparse_matrix),intent(in)    :: tAp,tAv,tApp
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(QU,2)) :: QP_current,QU_current
    integer                   :: i

    QP(n_time_step,:)=0.0
    QU(n_time_step,:)=0.0
    DP(n_time_step,:)=0.0
!    DU_half(n_time_step,:)=0.0
    
    do i=n_time_step-1,0-1
       ! print*,'i',i

       QP_current=QP(i+1,:)
       QU_current=QU(i+1,:)
       call LF_forward(QP_current,QU_current,tAv,tAp,tApp,DP(i+1,:),DU_half(i+1,:))
       QP(i,:)=QP_current
       QU(i,:)=QU_current
    end do
  end subroutine test_LF_backward

    subroutine test_RK4_forward(P,U,FP,FU,FP_half,FU_half,Ap,Av,App,n_time_step,&
                                size_v,nb_elem,DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: U
    real,dimension(0:n_time_step,size_v),intent(inout) :: P
    real,dimension(0:n_time_step,size_v),intent(inout) :: FP
    real,dimension(0:n_time_step,size_v),intent(inout) :: FU
    real,dimension(1:n_time_step,size_v),intent(inout) :: FP_half
    real,dimension(1:n_time_step,size_v),intent(inout) :: FU_half
    type(sparse_matrix),intent(in)    :: Ap,Av,App
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(U,2)) :: P_current,U_current
    integer                   :: i
    
    P(0,:)=0.0
    U(0,:)=0.0
    FP(0,:)=0.0
    FU(0,:)=0.0

    do i=1,n_time_step
       ! print*,'i',i
       
       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call RK4_forward(P_current,U_current,Ap,Av,App,                          &
            FP(i-1,:),FP_half(i,:),   &
            FP(i,:),FU(i-1,:),        &
            FU_half(i,:),FU(i,:),     &
            0.0*FP(i,:),0.0*FU(i,:))


       ! call RK4_forward(P_current,U_current,Ap,Av,App,                          &
       !      0.0*FP(i-1,:),0.0*FP_half(i,:),   &
       !      0.0*FP(i,:),0.0*FU(i-1,:),        &
       !      0.0*FU_half(i,:),0.0*FU(i,:),     &
       !      FP(i,:),FU(i,:))

       
       P(i,:)=P_current
       U(i,:)=U_current

    end do
  end subroutine test_RK4_forward

  subroutine test_RK4_forward2(P,U,FP,FU,FP_half,FU_half,Ap,Av,App,n_time_step,&
                                size_v,nb_elem,DoF,dx,GP,GU)
    real,dimension(0:n_time_step,size_v),intent(inout) :: U
    real,dimension(0:n_time_step,size_v),intent(inout) :: P
    real,dimension(0:n_time_step,size_v),intent(inout) :: FP
    real,dimension(0:n_time_step,size_v),intent(inout) :: FU
    real,dimension(0:n_time_step,size_v),intent(inout) :: GP
    real,dimension(0:n_time_step,size_v),intent(inout) :: GU
    real,dimension(1:n_time_step,size_v),intent(inout) :: FP_half
    real,dimension(1:n_time_step,size_v),intent(inout) :: FU_half
    type(sparse_matrix),intent(in)    :: Ap,Av,App
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(U,2)) :: P_current,U_current
    integer                   :: i
    
    P(0,:)=0.0
    U(0,:)=0.0
    FP(0,:)=0.0
    FU(0,:)=0.0

    do i=1,n_time_step
       ! print*,'i',i
       
       P_current=P(i-1,:)
       U_current=U(i-1,:)
       ! call RK4_forward(P_current,U_current,Ap,Av,App,                          &
       !      FP(i-1,:),FP_half(i,:),   &
       !      FP(i,:),FU(i-1,:),        &
       !      FU_half(i,:),FU(i,:),     &
       !      0.0*FP(i,:),0.0*FU(i,:))


       call RK4_forward(P_current,U_current,Ap,Av,App,                          &
            0.0*FP(i-1,:),0.0*FP_half(i,:),   &
            0.0*FP(i,:),0.0*FU(i-1,:),        &
            0.0*FU_half(i,:),0.0*FU(i,:),     &
            GP(i,:),GU(i,:))

       
       P(i,:)=P_current
       U(i,:)=U_current

    end do
 
  end subroutine test_RK4_forward2

  subroutine test_RK4_backward(QP,QU,DP,DU,tAv,tAp,tApp, &
                               n_time_step,size_v,nb_elem,DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: QU
    real,dimension(0:n_time_step,size_v),intent(inout) :: QP
    real,dimension(0:n_time_step,size_v),intent(inout) :: DP
    real,dimension(0:n_time_step,size_v),intent(inout) :: DU
    type(sparse_matrix),intent(in)    :: tAp,tAv,tApp
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx
    real,dimension(size_v) :: zero
    real,dimension(size(QP,2)) :: QP_current,QU_current
    integer                    :: i


    zero=0.0
    QP(n_time_step,:)=0.0
    QU(n_time_step,:)=0.0
    DP(n_time_step,:)=0.0
    !DP_half(n_time_step,:)=0.0
    DU(n_time_step,:)=0.0
    !DU_half(n_time_step,:)=0.0
    
    do i=n_time_step-1,0,-1
       ! print*,'i',i
       QP_current=QP(i+1,:)
       QU_current=QU(i+1,:)
       call RK4_forward(QP_current,QU_current,tAv,tAp,tApp,zero,zero,zero,zero,zero,zero,DP(i,:),DU(i,:))
!       call RK4_forward(QP_current,QU_current,tAv,tAp,tApp,zero,zero,zero,zero,zero,zero,DP(i+1,:),DU(i+1,:))
       QP(i,:)=QP_current
       QU(i,:)=QU_current
    end do
    
  end subroutine test_RK4_backward

  

  subroutine test_AB3_forward(P,U,FP,FU,Ap,Av,App,n_time_step,size_v,nb_elem,   &
                              DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: U
    real,dimension(0:n_time_step,size_v),intent(inout) :: P
    real,dimension(0:n_time_step,size_v),intent(inout)    :: FP
    real,dimension(0:n_time_step,size_v),intent(inout)    :: FU
    
    type(sparse_matrix),intent(in)    :: Ap,Av,App
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(U,2)) :: P_current,U_current
    real,dimension(size(U,2)) :: Pk1,Pk2
    real,dimension(size(U,2)) :: Uk1,Uk2
    integer                   :: i


    P(0:1,:)=0.0
    U(0:1,:)=0.0
    FP(0:1,:)=0.0
    FU(0:1,:)=0.0
    
    Pk1=0.0
    Pk2=0.0
    Uk1=0.0
    Uk2=0.0
    
    do i=2,n_time_step
       ! print*,'i',i

       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call AB3_forward(P_current,U_current,Ap,Av,App,FP(i,:),FU(i,:),      &
                        Pk1,Pk2,Uk1,Uk2)
       P(i,:)=P_current
       U(i,:)=U_current
    end do
  end subroutine test_AB3_forward

  subroutine test_AB3_backward(QP,QU,DP,DU,tAv,tAp,tApp,n_time_step,size_v,     &
                               nb_elem,DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: QU
    real,dimension(0:n_time_step,size_v),intent(inout) :: QP
    real,dimension(0:n_time_step,size_v),intent(inout) :: DP
    real,dimension(0:n_time_step,size_v),intent(inout) :: DU
    
    type(sparse_matrix),intent(in)    :: tAp,tAv,tApp
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(QU,2)) :: QP_current,QU_current
    real,dimension(size(QU,2)) :: QPk1,QPk2
    real,dimension(size(QU,2)) :: QUk1,QUk2
    integer                    :: i

    QPk1=0.0
    QPk2=0.0
    QUk1=0.0
    QUk2=0.0

   QP(n_time_step-1:n_time_step,:)=0.0
   QU(n_time_step-1:n_time_step,:)=0.0
   DP(n_time_step-1:n_time_step,:)=0.0
   DU(n_time_step-1:n_time_step,:)=0.0
     
   do i=n_time_step-2,0,-1
      
       QP_current=QP(i+1,:)
       QU_current=QU(i+1,:)
       call AB3_forward(QP_current,QU_current,tAv,tAp,tApp,DP(i+1,:),DU(i+1,:), &
                        QPk1,QPk2,QUk1,QUk2)
       QP(i,:)=QP_current
       QU(i,:)=QU_current
    end do
  end subroutine test_AB3_backward

  function inner_product(P,U,DP,DU)
    real,dimension(:,:),intent(in) :: U,P
    real,dimension(:,:),intent(in) :: DU,DP
    real                           :: inner_product
    integer                        :: i,j
    integer                        :: n_time_step

    n_time_step=size(U,1)
    inner_product=0.0

    
    do i=1,n_time_step
       do j=1,size(U,2)
          inner_product=inner_product+P(i,j)*DP(i,j)+U(i,j)*DU(i,j)
       end do
    end do
  end function inner_product
  
  function inner_product1(P,U,DP,DU)
    real,dimension(:,:),intent(in) :: U,P
    real,dimension(:,:),intent(in) :: DU,DP
    real                           :: inner_product1
    integer                        :: i,j
    integer                        :: n_time_step
    real,dimension(size(U,2)) :: EP,EU
    real :: b0,b1,b2,b3

    b0=55.0/24.0
    b1=-59.0/24.0
    b2=37.0/24.0
    b3=-9.0/24.0
    
    n_time_step=size(U,1)
    inner_product1=0.0

    
    do i=4,n_time_step-4

       
    EP=b0*DP(n_time_step-i+1,:)    ! &
        ! +b1*DP(n_time_step-i+2,:)  &
        ! +b2*DP(n_time_step-i+3,:)  &
        ! +b3*DP(n_time_step-i+4,:)

    EU=b0*DU(n_time_step-i+1,:)    ! &
        ! +b1*DU(n_time_step-i+2,:)  &
        ! +b2*DU(n_time_step-i+3,:)  &
        ! +b3*DU(n_time_step-i+4,:)


       
       do j=1,size(U,2)
          inner_product1=inner_product1+U(i,j)*EU(j)+P(i,j)*EP(j)
       end do
    end do
  end function inner_product1

  function inner_product2(FP,FU,QP,QU)
    real,dimension(:,:),intent(in) :: FU,FP
    real,dimension(:,:),intent(in) :: QU,QP
    real                           :: inner_product2
    integer                        :: i,j
    integer                        :: n_time_step
    real,dimension(size(FU,2)) :: EP,EU
    real :: b0,b1,b2,b3

    b0=55.0/24.0
    b1=-59.0/24.0
    b2=37.0/24.0
    b3=-9.0/24.0
    
    n_time_step=size(FU,1)
    inner_product2=0.0

    
    do i=4,n_time_step-4

       
    EP=b0*FP(i,:)       &
         +b1*FP(i-1,:)  &
         +b2*FP(i-2,:)  &
         +b3*FP(i-3,:)

    EU=b0*FU(i,:)       &
         +b1*FU(i-1,:)  &
         +b2*FU(i-2,:)  &
         +b3*FU(i-3,:)


       
       do j=1,size(FU,2)
          inner_product2=inner_product2+EU(j)*QU(n_time_step-i+1,j) &
                                     +EP(j)*QP(n_time_step-i+1,j)
       end do
    end do
  end function inner_product2
  

end module m_adjoint_test
