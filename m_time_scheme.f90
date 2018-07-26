module m_time_scheme
  use m_matrix
  implicit none

  contains

  subroutine LF_forward(P,U,Ap,Av,App,FP,FU)
    real,dimension(:)  ,intent(inout) :: P
    real,dimension(:)  ,intent(inout) :: U
    type(sparse_matrix),intent(in)    :: Ap
    type(sparse_matrix),intent(in)    :: Av
    type(sparse_matrix),intent(in)    :: App
    real,dimension(:)  ,intent(in)    :: FP
    real,dimension(:)  ,intent(in)    :: FU

    U=U-sparse_matmul(Av,P)-FU
    P=-sparse_matmul(App,P)-sparse_matmul(Ap,U)-FP
  end subroutine LF_forward
  

  subroutine RK4_forward(P,U,Ap,Av,App,FP_0,FP_half,FP_1,FU_0,FU_half,FU_1)
    real,dimension(:)  ,intent(inout) :: P
    real,dimension(:)  ,intent(inout) :: U
    type(sparse_matrix),intent(in)    :: Ap
    type(sparse_matrix),intent(in)    :: Av
    type(sparse_matrix),intent(in)    :: App
    real,dimension(:)  ,intent(in)    :: FP_0,FP_half,FP_1
    real,dimension(:)  ,intent(in)    :: FU_0,FU_half,FU_1

    real,dimension(size(P)) :: Uk1,Uk2,Uk3,Uk4
    real,dimension(size(P)) :: Pk1,Pk2,Pk3,Pk4

    Uk1=-sparse_matmul(Av,P)-FU_0
    Pk1=-sparse_matmul(App,P)-sparse_matmul(Ap,U)-FP_0
    
    Uk2=-sparse_matmul(Av,P+0.5*Pk1)-FU_half
    Pk2=-sparse_matmul(App,P+0.5*Pk1)-sparse_matmul(Ap,U+0.5*Uk1)-FP_half

    Uk3=-sparse_matmul(Av,P+0.5*Pk2)-FU_half
    Pk3=-sparse_matmul(App,P+0.5*Pk2)-sparse_matmul(Ap,U+0.5*Uk2)-FP_half

    Uk4=-sparse_matmul(Av,P+Pk3)-FU_1
    Pk4=-sparse_matmul(App,P+Pk3)-sparse_matmul(Ap,U+Uk3)-FP_1

    P=P+(1.0/6.0)*(Pk1+2.0*Pk2+2.0*Pk3+Pk4)
    U=U+(1.0/6.0)*(Uk1+2.0*Uk2+2.0*Uk3+Uk4)
  end subroutine RK4_forward


  subroutine AB3_forward(P,U,Ap,Av,App,FP,FU,Pk1,Pk2,Pk3,Uk1,Uk2,Uk3)
    real,dimension(:)  ,intent(inout) :: P
    real,dimension(:)  ,intent(inout) :: U
    type(sparse_matrix),intent(in)    :: Ap
    type(sparse_matrix),intent(in)    :: Av
    type(sparse_matrix),intent(in)    :: App
    real,dimension(:)  ,intent(in)    :: FP,FU
    real,dimension(:)  ,intent(inout) :: Pk1,Pk2,Pk3
    real,dimension(:)  ,intent(inout) :: Uk1,Uk2,Uk3
    
    real,dimension(size(P)) :: Pk0,Uk0

    Uk0=-sparse_matmul(Av,P)-FU
    Pk0=-sparse_matmul(App,P)-sparse_matmul(Ap,U)-FP

    P=P+55.0/24.0*Pk0         &
       -59.0/24.0*Pk1         &
       +37.0/24.0*Pk2         &
         -9.0/24.0*Pk3

    U=U+55.0/24.0*Uk0         &
       -59.0/24.0*Uk1         &
       +37.0/24.0*Uk2         &
       -9.0/24.0*Uk3

    Pk3=Pk2
    Uk3=Uk2
    Pk2=Pk1
    Uk2=Uk1
    Pk1=Pk0
    Uk1=Uk0
    end subroutine AB3_forward
  
end module m_time_scheme
