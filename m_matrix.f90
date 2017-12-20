module m_matrix
  implicit none
  
  type sparse_matrix
     integer :: NNN
     integer :: nb_ligne
     real,dimension(:),allocatable :: Values
     integer,dimension(:),allocatable :: IA
     integer,dimension(:),allocatable :: JA
  end type sparse_matrix

contains

  function det2(A)
    real,dimension(2,2),intent(in) :: A
    real                           :: det2

    det2=A(1,1)*A(2,2)-A(1,2)*A(2,1)
  end function det2


  function det3(A)
    real,dimension(3,3),intent(in) :: A
    real                           :: det3

    det3=A(1,1)*det2(A(2:3,2:3))           &
    -A(2,1)*(A(1,2)*A(3,3)-A(1,3)*A(3,2))  &
    +A(3,1)*det2(A(1:2,2:3))
    
  end function det3


  function det4(A)
    real,dimension(4,4),intent(in) :: A
    real                           :: det4
    
    det4=A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+&
         A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+&
         A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))  &
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+&
         A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
         A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+&
         A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+&
         A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+&
         A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+&
         A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
  end function det4


  function det5(A)
    real,dimension(5,5),intent(in) :: A
    real                           :: det5
    integer                        :: i
    real,dimension(4,4)            :: B

    do i=1,5

       B(1:i-1,:)=A(1:i-1,2:5)
       B(i:4,:)=A(i+1:5,2:5)

       det5=det5+(-1)**(i+1)*A(i,1)*det4(B)

    end do

  end function det5

  function inv5(A)
    real,dimension(5,5),intent(in) :: A
    real,dimension(5,5)            :: inv5
    real,dimension(4,4)            :: B
    integer                        :: i,i1,j

    inv5=0.0
    B=0.0
    
    do i=1,5
       do j=1,5
          B(1:i-1,1:j-1)=A(1:i-1,1:j-1)
          B(i:4,1:j-1)=A(i+1:5,1:j-1)
          B(1:i-1,j:4)=A(1:i-1,j+1:5)
          B(i:4,j:4)=A(i+1:5,j+1:5)

          inv5(i,j)=(-1)**(i+j)*det4(B)

       end do
    end do
    inv5=transpose(inv5)
    inv5=(1.0/det5(A))*inv5
  end function inv5

  subroutine init_sparse_matrix(A,NNN,nb_ligne)
    type(sparse_matrix),intent(inout) :: A
    integer            ,intent(in)    :: NNN
    integer            ,intent(in)    :: nb_ligne

    A%NNN=NNN
    A%nb_ligne=nb_ligne
    allocate(A%Values(NNN))
    allocate(A%IA(0:nb_ligne))
    allocate(A%JA(1:NNN))
    
  end subroutine init_sparse_matrix

  function get_NNN(A)
    real,dimension(:,:),intent(in) :: A
    integer                        :: get_NNN
    integer                        :: i,j

    get_NNN=0
    
    do i=1,size(A,1)
       do j=1,size(A,2)

          if(A(i,j).ne.0.0) then
             get_NNN=get_NNN+1
          end if
       end do
    end do
  end function get_NNN

  subroutine Full2Sparse(Full,Sparse)
    real,dimension(:,:),intent(in)  :: Full
    type(sparse_matrix),intent(out) :: Sparse

    integer :: i,j,ivalues,jja,NNN_ligne
    
    call init_sparse_matrix(Sparse,get_NNN(Full),size(Full,1))
    
    ivalues=1
    jja=1

    do i=1,size(Full,1)
       NNN_ligne=0
       do j=1,size(Full,2)

          if (Full(i,j).ne.0) then
             Sparse%Values(ivalues)=Full(i,j)
             ivalues=ivalues+1
             NNN_ligne=NNN_ligne+1
             Sparse%JA(jja)=j
             jja=jja+1
          end if
       end do
       Sparse%IA(i)=Sparse%IA(i-1)+NNN_ligne
    end do

  end subroutine Full2Sparse

  subroutine print_sparse_matrix(A)
    type(sparse_matrix),intent(in) :: A

    print*,'NNN :',A%NNN
    print*,'nb_ligne :',A%nb_ligne
    print*,'Values :',A%Values
    print*,'IA :',A%IA
    print*,'JA :',A%JA
    
  end subroutine print_sparse_matrix

  function sparse_matmul(A,X)
    type(sparse_matrix),intent(in) :: A
    real,dimension(:)  ,intent(in) :: X
    real,dimension(size(X))        :: sparse_matmul

    integer :: i,j


    print*,'test',size(X),size(A%Values),A%nb_ligne
    call print_sparse_matrix(A)
    sparse_matmul=0.0

    do i=1,A%nb_ligne
       do j=A%IA(i-1)+1,A%IA(i)
          sparse_matmul(i)=sparse_matmul(i)+A%Values(j)*X(A%JA(j))
       end do
    end do
  end function sparse_matmul

  
end module m_matrix
