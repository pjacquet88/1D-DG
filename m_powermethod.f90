! This module contain the power method used to determine the maximal eigen values
! employed to determine the CFL
module m_powermethod
  use m_matrix

  implicit none

  public  :: power_method,power_method_sparse
  private :: norm 

contains

  !*************** PUBLIC *******************************************************
  
  ! Power method employed on a full matrix
  subroutine power_method(A,max_value,k_max,epsilon)
    real,dimension(:,:),intent(in)  :: A
    real               ,intent(out) :: max_value
    integer            ,intent(in)  :: k_max
    real               ,intent(in)  :: epsilon
    real,dimension(:),allocatable   :: z,q1,q2
    real                            :: r
    integer                         :: k
    integer :: i

    r=1.0
    k=1

    allocate(z(size(A,1)),q1(size(A,1)),q2(size(A,1)))

    do i=1,size(q1)
       q1(i)=rand()
    end do
    q1=1/(norm(q1))*q1
    q2=q1

    do while((r.ge.epsilon).and.(k.le.k_max))
       z=matmul(A,q2)
       q2=1.0/(norm(z))*z
       max_value=dot_product(q2,matmul(A,q2))
       r=norm(q2-q1)
       q1=q2
       k=k+1
    end do
    deallocate(z,q1,q2)
  end subroutine power_method


  ! Power method employed on asparse matrix
  subroutine power_method_sparse(A,max_value,k_max,epsilon)
    type(sparse_matrix),intent(in)  :: A
    real               ,intent(out) :: max_value
    integer            ,intent(in)  :: k_max
    real               ,intent(in)  :: epsilon
    real,dimension(:),allocatable   :: z,q1,q2
    real                            :: r
    integer                         :: k
    r=1.0
    k=1

    allocate(z(A%nb_ligne),q1(A%nb_ligne),q2(A%nb_ligne))
    q1=1.0
    q2=1.0
    do while((r.ge.epsilon).and.(k.le.k_max))
       z=sparse_matmul(A,q2)
       q2=1.0/(norm(z))*z
       max_value=dot_product(q2,sparse_matmul(A,q2))
       r=norm(q2-q1)/norm(q1)
       q1=q2
       k=k+1
    end do
    deallocate(z,q1,q2)
 end subroutine power_method_sparse


 !********** PRIVATE ************************************************************
 
 ! Calculate the norm of a vector
 function norm(v)
    real,dimension(:),intent(in) :: v
    real                         :: norm
    norm=sqrt(dot_product(v,v))
  end function norm

end module m_powermethod
