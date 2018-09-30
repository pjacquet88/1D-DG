module m_matrix

  implicit none
  
  type sparse_matrix
     integer :: NNN
     integer :: nb_ligne
     real,dimension(:),allocatable :: Values
     integer,dimension(:),allocatable :: IA
     integer,dimension(:),allocatable :: JA
  end type sparse_matrix


  public  :: sparse_matrix,                                                     &
             LU_inv,                                                            &
             init_sparse_matrix,get_NNN,Full2Sparse,print_sparse_matrix,        &
             sparse_matmul,is_sym,free_sparse_matrix,sparse_MV,sparse_MM

    interface sparse_matmul
     module procedure sparse_MV
     module procedure sparse_MM
  end interface sparse_matmul


contains


  !************ Matrices Inversion *******************************
  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  function LU_inv(A)
    real,dimension(:,:), intent(in) :: A
    real,dimension(size(A,1),size(A,2)) :: LU_inv
  
    real,   dimension(size(A,1)) :: work  ! work array for LAPACK
    integer,dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info,i

    ! External procedures defined in LAPACK
    external SGETRF
    external SGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    LU_inv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.+

    ! call SGETRF(n, n, LU_inv, n, ipiv, info)
    call DGETRF(n, n, LU_inv, n, ipiv, info)

    if (info /= 0) then
       ! do i=1,n
       !    write(22,*) A(i,:)
       ! end do
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    ! call SGETRI(n, LU_inv, n, ipiv, work, n, info)
    call DGETRI(n, LU_inv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if

  end function LU_inv


  !************ UNCTION SPARSE MATRIX ******************************************
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

  
  subroutine free_sparse_matrix(A)
       type(sparse_matrix),intent(inout) :: A
    deallocate(A%Values)
    deallocate(A%IA)
    deallocate(A%JA)
  end subroutine free_sparse_matrix
    

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
    integer                         :: i,j,ivalues,jja,NNN_ligne
    
    call init_sparse_matrix(Sparse,get_NNN(Full),size(Full,1))
    ivalues=1
    jja=1

    do i=1,size(Full,1)
       NNN_ligne=0
       do j=1,size(Full,2)

          if (Full(i,j).ne.0.0) then
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

  subroutine Sparse2Full(sparse,full)
    type(sparse_matrix),intent(in)    :: sparse
    real,dimension(:,:),allocatable,intent(out) :: full

    integer :: i,j,nb_value,iia
    
    allocate(full(sparse%nb_ligne,sparse%nb_ligne))
    full=0.0

    do i=1,sparse%nb_ligne
       nb_value=sparse%IA(i-1)
       do iia=1,sparse%IA(i)-sparse%IA(i-1)
          full(i,sparse%JA(nb_value+iia))=sparse%Values(nb_value+iia)
       end do
    end do
  end subroutine Sparse2Full
    
  subroutine transpose_sparse(sparse,t_sparse)
    type(sparse_matrix),intent(in)  :: sparse
    type(sparse_matrix),intent(out) :: t_sparse
    real,dimension(:,:),allocatable :: full

    call init_sparse_matrix(t_sparse,sparse%NNN,sparse%nb_ligne)

    call Sparse2Full(sparse,full)
    call Full2Sparse(transpose(full),t_sparse)
  end subroutine transpose_sparse


  function sparse_MV(A,X)
    type(sparse_matrix),intent(in) :: A
    real,dimension(:)  ,intent(in) :: X
    real,dimension(size(X))        :: sparse_MV
    integer                        :: i,j

    sparse_MV=0.0
    do i=1,A%nb_ligne
       do j=A%IA(i-1)+1,A%IA(i)
          sparse_MV(i)=sparse_MV(i)+A%Values(j)*X(A%JA(j))
       end do
    end do
  end function sparse_MV

  function sparse_MM(A,B)
    type(sparse_matrix),intent(in)  :: A
    type(sparse_matrix),intent(in)  :: B
    type(sparse_matrix)             :: sparse_MM
    real,dimension(:,:),allocatable :: A_full,B_full,C_full

    call Sparse2Full(A,A_full)
    call Sparse2Full(B,B_full)
    C_full=matmul(A_full,B_full)
    call Full2Sparse(C_full,sparse_MM)

  end function sparse_MM


  !************ TEST FUNCTIONS  *************************************************

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
