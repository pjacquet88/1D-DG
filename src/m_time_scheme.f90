! This module contains the differnet time scheme used

module m_time_scheme
  use m_matrix
  implicit none

  public :: LF_forward, RK4_forward, AB3_forward

  contains

  ! Increments P and U using one Leap-Frog time step
  subroutine LF_forward(P,U,Ap,Av,App,FP,FU)
    real(mp),dimension(:),intent(inout) :: P    ! pressure
    real(mp),dimension(:),intent(inout) :: U    ! velocity
    type(sparse_matrix)  ,intent(in)    :: Ap
    type(sparse_matrix)  ,intent(in)    :: Av
    type(sparse_matrix)  ,intent(in)    :: App
    real(mp),dimension(:),intent(in)    :: FP  ! RHS of the dtP equation
    real(mp),dimension(:),intent(in)    :: FU  ! RHS of the dtV equation

    U=U-sparse_matmul(Av,P)-FU
    P=-sparse_matmul(App,P)-sparse_matmul(Ap,U)-FP
  end subroutine LF_forward


  ! Increments P and U using one RK4 time step
  subroutine RK4_forward(P,U,Ap,Av,App,FP_0,FP_half,FP_1,FU_0,FU_half,FU_1,GP,GU)
    real(mp),dimension(:),intent(inout) :: P         ! pressure
    real(mp),dimension(:),intent(inout) :: U         ! velocity
    type(sparse_matrix)  ,intent(in)    :: Ap
    type(sparse_matrix)  ,intent(in)    :: Av
    type(sparse_matrix)  ,intent(in)    :: App
    real(mp),dimension(:),intent(in)    :: FP_0     ! dtp equation RHS at t=n
    real(mp),dimension(:),intent(in)    :: FP_half  ! dtp equation RHS at t=n+1/2
    real(mp),dimension(:),intent(in)    :: FP_1     ! dtp equation RHS at t=n+1
    real(mp),dimension(:),intent(in)    :: FU_0     ! dtu equation RHS at t=n
    real(mp),dimension(:),intent(in)    :: FU_half  ! dtu equation RHS at t=n+1/2
    real(mp),dimension(:),intent(in)    :: FU_1     ! dtu equation RHS at t=n+1
    real(mp),dimension(:),optional      :: GP,GU    ! optional : Adjoint RHS

    real(mp),dimension(size(P)) :: Uk1,Uk2,Uk3,Uk4
    real(mp),dimension(size(P)) :: Pk1,Pk2,Pk3,Pk4

    Uk1=-sparse_matmul(Av,P)+FU_0
    Pk1=-sparse_matmul(App,P)-sparse_matmul(Ap,U)+FP_0

    Uk2=-sparse_matmul(Av,P+0.5*Pk1)+FU_half
    Pk2=-sparse_matmul(App,P+0.5*Pk1)-sparse_matmul(Ap,U+0.5*Uk1)+FP_half

    Uk3=-sparse_matmul(Av,P+0.5*Pk2)+FU_half
    Pk3=-sparse_matmul(App,P+0.5*Pk2)-sparse_matmul(Ap,U+0.5*Uk2)+FP_half

    Uk4=-sparse_matmul(Av,P+Pk3)+FU_1
    Pk4=-sparse_matmul(App,P+Pk3)-sparse_matmul(Ap,U+Uk3)+FP_1

    if ((present(GP)).and.(present(GU))) then
       P=P+(1.0/6.0)*(Pk1+2.0*Pk2+2.0*Pk3+Pk4)+GP
       U=U+(1.0/6.0)*(Uk1+2.0*Uk2+2.0*Uk3+Uk4)+GU
    else
       P=P+(1.0/6.0)*(Pk1+2.0*Pk2+2.0*Pk3+Pk4)
       U=U+(1.0/6.0)*(Uk1+2.0*Uk2+2.0*Uk3+Uk4)
    end if
  end subroutine RK4_forward


  ! Increments P and U using a AB3 time step
  subroutine AB3_forward(P,U,Ap,Av,App,FP0,FU0,FP1,FU1,FP2,FU2,                 &
                         Pk1,Pk2,Uk1,Uk2,GP,GU)
    real(mp),dimension(:),intent(inout) :: P         ! pressure
    real(mp),dimension(:),intent(inout) :: U         ! velocity
    type(sparse_matrix)  ,intent(in)    :: Ap
    type(sparse_matrix)  ,intent(in)    :: Av
    type(sparse_matrix)  ,intent(in)    :: App
    real(mp),dimension(:),intent(in)    :: FP0,FU0   ! dtP and dtU RHS at t=n
    real(mp),dimension(:),intent(in)    :: FP1,FU1   ! dtP and dtU RHS at t=n-1
    real(mp),dimension(:),intent(in)    :: FP2,FU2   ! dtP and dtU RHS at t=n-2
    real(mp),dimension(:),intent(inout) :: Pk1,Pk2   ! dtP evaluates at t=n-1/n-2
    real(mp),dimension(:),intent(inout) :: Uk1,Uk2   ! dtU evaluates at t=n-1/n-2
    real(mp),dimension(:),optional      :: GP,GU     ! Adjoint RHS

    real(mp),dimension(size(P)) :: Pk0
    real(mp),dimension(size(U)) :: Uk0
    real(mp) :: b0,b1,b2

    if ((present(GP)).and.(present(GU))) then
       Uk0=-sparse_matmul(Av,P)
       Pk0=-sparse_matmul(App,P)-sparse_matmul(Ap,U)
    else
       Uk0=-sparse_matmul(Av,P)+FU0
       Pk0=-sparse_matmul(App,P)-sparse_matmul(Ap,U)+FP0
    end if

    P=P+23.0/12.0*Pk0         &
       -16.0/12.0*Pk1         &
       +5.0/12.0*Pk2

    U=U+23.0/12.0*Uk0         &
       -16.0/12.0*Uk1         &
       +5.0/12.0*Uk2

    Pk2=Pk1
    Uk2=Uk1
    Pk1=Pk0
    Uk1=Uk0

    if ((present(GP)).and.(present(GU))) then
       P=P+GP
       U=U+GU
    end if
  end subroutine AB3_forward

end module m_time_scheme
