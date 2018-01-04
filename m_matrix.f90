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

  function inv(A)
    real,dimension(:,:),intent(in)          :: A
    real,dimension(size(A,1),size(A,1))     :: inv
    real,dimension(size(A,1)-1,size(A,1)-1) :: B
    integer                                 :: i,i1,j,k

    inv=0.0
    B=0.0

    ! do k=1,size(A,1)
    !    print*,A(k,:)
    ! end do

    do i=1,size(A,1)
       do j=1,size(A,1)
          B(1:i-1,1:j-1)=A(1:i-1,1:j-1)
          B(i:size(A,1)-1,1:j-1)=A(i+1:size(A,1),1:j-1)
          B(1:i-1,j:size(A,1)-1)=A(1:i-1,j+1:size(A,1))
          B(i:size(A,1)-1,j:size(A,1)-1)=A(i+1:size(A,1),j+1:size(A,1))
          inv(i,j)=(-1)**(i+j)*FindDet(B)

          ! print*,'inversion B',i,j
          ! do k=1,size(B,1)
          !    print*,B(k,:)
          ! end do
          

       end do
    end do
    inv=transpose(inv)
    inv=(1.0/FindDet(A))*inv
  end function inv


  real function FindDet(A)
    implicit none
    real,intent(in),dimension(:,:) :: A
    real,dimension(size(A,1),size(A,2)) :: matrix
    integer              :: n
    real :: m, temp
    integer :: i, j, k, l
    logical :: detexists = .TRUE.

    matrix=A

    n=size(matrix,1)
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
       IF (matrix(k,k) == 0) THEN
          DetExists = .FALSE.
          DO i = k+1, n
             IF (matrix(i,k) /= 0) THEN
                DO j = 1, n
                   temp = matrix(i,j)
                   matrix(i,j)= matrix(k,j)
                   matrix(k,j) = temp
                END DO
                DetExists = .TRUE.
                l=-l
                EXIT
             ENDIF
          END DO
          IF (DetExists .EQV. .FALSE.) THEN
             FindDet = 0
             return
          END IF
       ENDIF
       DO j = k+1, n
          m = matrix(j,k)/matrix(k,k)
          DO i = k+1, n
             matrix(j,i) = matrix(j,i) - m*matrix(k,i)
          END DO
       END DO
    END DO

    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
       FindDet = FindDet * matrix(i,i)
    END DO

  END FUNCTION FindDet

  subroutine init_sparse_matrix(A,NNN,nb_ligne)
    type(sparse_matrix),intent(inout) :: A
    integer            ,intent(in)    :: NNN
    integer            ,intent(in)    :: nb_ligne

    A%NNN=NNN
    A%nb_ligne=nb_ligne
    allocate(A%Values(NNN))
    allocate(A%IA(0:nb_ligne))
    A%IA(0)=0
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

    sparse_matmul=0.0

    do i=1,A%nb_ligne
       do j=A%IA(i-1)+1,A%IA(i)
          sparse_matmul(i)=sparse_matmul(i)+A%Values(j)*X(A%JA(j))
       end do
    end do
  end function sparse_matmul
  
end module m_matrix
