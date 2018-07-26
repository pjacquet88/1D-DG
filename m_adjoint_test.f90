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
    allocate(test%FP_half(0:n_time_step,test%size_v),test%FU_half(0:n_time_step,test%size_v))
    allocate(test%DP(0:n_time_step,test%size_v),test%DU(0:n_time_step,test%size_v))
    allocate(test%DP_half(0:n_time_step,test%size_v),test%DU_half(0:n_time_step,test%size_v))
    
    call random_number(test%P)
    call random_number(test%U)
    call random_number(test%QP)
    call random_number(test%QU)
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
    
  end subroutine init_adjoint_test

  
  subroutine forward_test(test)
    type(t_adjoint_test),intent(inout) :: test
    
   test%U(0,:)=0.0
   test%P(0,:)=0.0
   test%FU(0,:)=0.0
   test%FP(0,:)=0.0

    if (test%time_scheme.eq.'RK4') then
       call test_RK4_forward(test%P,test%U,test%FP,test%FU,test%FP_half,        &
                             test%FU_half,test%Ap,test%Av,test%App,             &
                             test%n_time_step,test%size_v,test%nb_elem,test%DoF,&
                             test%dx)
       
    else if (test%time_scheme.eq.'LF') then
       call test_LF_forward(test%P,test%U,test%FP,test%FU_half,test%Ap,test%Av, &
                            test%App,test%n_time_step,test%size_v,test%nb_elem, &
                            test%DoF,test%dx)

    else if (test%time_scheme.eq.'AB3') then
       call  test_AB3_forward(test%P,test%U,test%FP,test%FU,test%Ap,test%Av,    &
                              test%App,test%n_time_step,test%size_v,            &
                              test%nb_elem,test%DoF,test%dx)
       
    else
       print*,'not recongnized time scheme'
    end if

  end subroutine forward_test
     
  subroutine backward_test(test)
    type(t_adjoint_test),intent(inout) :: test

    test%QU(0,:)=0.0
    test%QP(0,:)=0.0
    test%DU(0,:)=0.0
    test%DP(0,:)=0.0

    if (test%time_scheme.eq.'RK4') then
       call test_RK4_forward(test%QP,test%QU,test%DP,test%DU,test%DP_half,        &
            test%DU_half,test%tAv,test%tAp,test%tApp,             &
            test%n_time_step,test%size_v,test%nb_elem,test%DoF,&
            test%dx)

    else if (test%time_scheme.eq.'LF') then
       call test_LF_forward(test%QP,test%QU,test%DP,test%DU_half,test%tAv,test%tAp, &
            test%tApp,test%n_time_step,test%size_v,test%nb_elem, &
            test%DoF,test%dx)

    else if (test%time_scheme.eq.'AB3') then
       call  test_AB3_forward(test%QP,test%QU,test%DP,test%DU,test%tAv,test%tAp,    &
            test%tApp,test%n_time_step,test%size_v,            &
            test%nb_elem,test%DoF,test%dx)

    else
       print*,'not recongnized time scheme'
    end if

  end subroutine backward_test

  subroutine test_LF_forward(P,U,FP,FU_half,Ap,Av,App,n_time_step,       &
       size_v,nb_elem,DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: P
    real,dimension(0:n_time_step,size_v),intent(inout) :: U
    real,dimension(0:n_time_step,size_v),intent(in)    :: FP
    real,dimension(0:n_time_step,size_v),intent(in)    :: FU_half
    type(sparse_matrix),intent(in)    :: Ap,Av,App
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(U,2)) :: P_current,U_current
    integer                   :: i

    call print_vect(U(0,:),nb_elem,DoF,dx,.true.,0,'FU')
    call print_vect(P(0,:),nb_elem,DoF,dx,.true.,0,'FP')
    do i=1,n_time_step
       ! print*,'i',i

       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call LF_forward(P_current,U_current,Ap,Av,App,FP(i-1,:),FU_half(i-1,:))
       P(i,:)=P_current
       U(i,:)=U_current

       call print_vect(U(i,:),nb_elem,DoF,dx,.true.,i,'FU')
       call print_vect(P(i,:),nb_elem,DoF,dx,.true.,i,'FP')
    end do
  end subroutine test_LF_forward

    subroutine test_RK4_forward(P,U,FP,FU,FP_half,FU_half,Ap,Av,App,n_time_step,       &
       size_v,nb_elem,DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: U
    real,dimension(0:n_time_step,size_v),intent(inout) :: P
    real,dimension(0:n_time_step,size_v),intent(in)    :: FP
    real,dimension(0:n_time_step,size_v),intent(in)    :: FU
    real,dimension(0:n_time_step,size_v),intent(in)    :: FP_half
    real,dimension(0:n_time_step,size_v),intent(in)    :: FU_half
    type(sparse_matrix),intent(in)    :: Ap,Av,App
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(U,2)) :: P_current,U_current
    integer                   :: i

    call print_vect(U(0,:),nb_elem,DoF,dx,.true.,0,'FU')
    call print_vect(P(0,:),nb_elem,DoF,dx,.true.,0,'FP')
    do i=1,n_time_step
       ! print*,'i',i

       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call RK4_forward(P_current,U_current,Ap,Av,App,FP(i-1,:),FP_half(i-1,:),   &
            FP(i,:),FU(i-1,:),FU_half(i-1,:),FU(i,:))
       P(i,:)=P_current
       U(i,:)=U_current

       call print_vect(U(i,:),nb_elem,DoF,dx,.true.,i,'FU')
       call print_vect(P(i,:),nb_elem,DoF,dx,.true.,i,'FP')
    end do
  end subroutine test_RK4_forward

  subroutine test_AB3_forward(P,U,FP,FU,Ap,Av,App,n_time_step,size_v,nb_elem,   &
                              DoF,dx)
    real,dimension(0:n_time_step,size_v),intent(inout) :: U
    real,dimension(0:n_time_step,size_v),intent(inout) :: P
    real,dimension(0:n_time_step,size_v),intent(in)    :: FP
    real,dimension(0:n_time_step,size_v),intent(in)    :: FU
    
    type(sparse_matrix),intent(in)    :: Ap,Av,App
    integer            ,intent(in)    :: n_time_step
    integer            ,intent(in)    :: size_v
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    real               ,intent(in)    :: dx

    real,dimension(size(U,2)) :: P_current,U_current
    real,dimension(size(U,2)) :: Pk1,Pk2,Pk3
    real,dimension(size(U,2)) :: Uk1,Uk2,Uk3
    integer                   :: i

    call print_vect(U(0,:),nb_elem,DoF,dx,.true.,0,'FU')
    call print_vect(P(0,:),nb_elem,DoF,dx,.true.,0,'FP')

    Pk1=0.0
    Pk2=0.0
    Pk3=0.0
    Uk1=0.0
    Uk2=0.0
    Uk3=0.0
    
    do i=1,n_time_step
       ! print*,'i',i

       P_current=P(i-1,:)
       U_current=U(i-1,:)
       call AB3_forward(P_current,U_current,Ap,Av,App,FP(i-1,:),FU(i-1,:),      &
                        Pk1,Pk2,Pk3,Uk1,Uk2,Uk3)
       P(i,:)=P_current
       U(i,:)=U_current

       call print_vect(U(i,:),nb_elem,DoF,dx,.true.,i,'FU')
       call print_vect(P(i,:),nb_elem,DoF,dx,.true.,i,'FP')
    end do
  end subroutine test_AB3_forward

  
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
          inner_product=inner_product+U(i,j)*DU(n_time_step-i+1,j)                &
                                     +P(i,j)*DP(n_time_step-i+1,j)
       end do
    end do
  end function inner_product

end module m_adjoint_test
