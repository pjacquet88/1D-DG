module m_matrix

  implicit none
  
  type sparse_matrix
     integer :: NNN
     integer :: nb_ligne
     real,dimension(:),allocatable :: Values
     integer,dimension(:),allocatable :: IA
     integer,dimension(:),allocatable :: JA
  end type sparse_matrix

  public  :: sparse_matrix,                                                    &
             LU_inv,                                                           &
             init_sparse_matrix,get_NNN,Full2Sparse,print_sparse_matrix,       &
             sparse_matmul,is_sym

! private :: FindDet

contains


  !************ INVERSION DIRECTE PETITE MATRICE *******************************

  
  
  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  function LU_inv(A)
    real,dimension(:,:), intent(in) :: A
    real,dimension(size(A,1),size(A,2)) :: LU_inv
    real,dimension(size(A,1),size(A,2)) :: test

    real,   dimension(size(A,1)) :: work  ! work array for LAPACK
    integer,dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info,i

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    LU_inv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.+

    call SGETRF(n, n, LU_inv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call SGETRI(n, LU_inv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if

    test=matmul(A,LU_inv)
    
  end function LU_inv

  
  ! function inv(A)
  !   real,dimension(:,:),intent(in)          :: A
  !   real,dimension(size(A,1),size(A,1))     :: inv,test
  !   real,dimension(size(A,1)-1,size(A,1)-1) :: B
  !   integer                                 :: i,i1,j,k

  !   inv=0.0
  !   B=0.0

  !   do i=1,size(A,1)
  !      do j=1,size(A,1)
  !         B(1:i-1,1:j-1)=A(1:i-1,1:j-1)
  !         B(i:size(A,1)-1,1:j-1)=A(i+1:size(A,1),1:j-1)
  !         B(1:i-1,j:size(A,1)-1)=A(1:i-1,j+1:size(A,1))
  !         B(i:size(A,1)-1,j:size(A,1)-1)=A(i+1:size(A,1),j+1:size(A,1))
  !         inv(i,j)=(-1)**(i+j)*FindDet(B)
  !      end do
  !   end do
  !   inv=transpose(inv)
  !   inv=(1.0/FindDet(A))*inv
  !   print*,'Det B2L :',FindDet(A),'ordre :',size(A,1)-1

    
  !   test=matmul(A,inv)
    
  ! end function inv


  ! real function FindDet(A)
  !   implicit none
  !   real,intent(in),dimension(:,:) :: A
  !   real,dimension(size(A,1),size(A,2)) :: matrix
  !   integer              :: n
  !   real :: m, temp
  !   integer :: i, j, k, l
  !   logical :: detexists = .TRUE.

  !   matrix=A

  !   n=size(matrix,1)
  !   l = 1
  !   !Convert to upper triangular form
  !   DO k = 1, n-1
  !      IF (matrix(k,k) == 0) THEN
  !         DetExists = .FALSE.
  !         DO i = k+1, n
  !            IF (matrix(i,k) /= 0) THEN
  !               DO j = 1, n
  !                  temp = matrix(i,j)
  !                  matrix(i,j)= matrix(k,j)
  !                  matrix(k,j) = temp
  !               END DO
  !               DetExists = .TRUE.
  !               l=-l
  !               EXIT
  !            ENDIF
  !         END DO
  !         IF (DetExists .EQV. .FALSE.) THEN
  !            FindDet = 0
  !            return
  !         END IF
  !      ENDIF
  !      DO j = k+1, n
  !         m = matrix(j,k)/matrix(k,k)
  !         DO i = k+1, n
  !            matrix(j,i) = matrix(j,i) - m*matrix(k,i)
  !         END DO
  !      END DO
  !   END DO

  !   !Calculate determinant by finding product of diagonal elements
  !   FindDet = l
  !   DO i = 1, n
  !      FindDet = FindDet * matrix(i,i)
  !   END DO

  ! END FUNCTION FindDet


  !************FUNCTION SPARSE MATRIX ******************************************
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



  !************ FONCTIONS TESTS ************************************************

  subroutine print_sparse_matrix(A)
    type(sparse_matrix),intent(in) :: A

    print*,'NNN :',A%NNN
    print*,'nb_ligne :',A%nb_ligne
    print*,'Values :',A%Values
    print*,'IA :',A%IA
    print*,'JA :',A%JA
    
  end subroutine print_sparse_matrix

  
  subroutine is_sym(A)
  real,dimension(:,:),intent(in) :: A
  integer :: i,j
  logical :: test

  test=.TRUE.
  
  do i=1,size(A,1)
     do j=i,size(A,1)

        if (A(i,j).ne.A(j,i)) then
           test=.FALSE.
           exit
        end if
     end do
     if (.not.test) exit
  end do

  if (test) then
     print*,'IS SYMETRIC'
  else
     print*,'IS NOT SYMETRIC'
  end if
end subroutine is_sym
  
end module m_matrix
