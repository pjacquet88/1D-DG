module m_polynom
  use m_matrix
  implicit none

  real               ,dimension(:,:),allocatable :: B2L,L2B,D0,D1
  real               ,dimension(:,:),allocatable :: base_b
  real               ,dimension(:,:),allocatable :: base_l
  type(sparse_matrix)                            :: D0_sparse,D1_sparse
  

  public  :: eval_polynom_b,eval_polynom_l,                                     &
             init_basis_b,init_basis_l,                                         &
             create_L2B,create_B2L,create_derive,                               &
             free_L2B,free_B2L,free_derive,                                     &
             Bernstein2Lagrange,Lagrange2Bernstein,B2L,L2B,D1,D0
  
  private :: D0_sparse,D1_sparse,                                               &
       C,b_basis,l_basis
contains


  !****************FONCTIONS ANNEXES**********************************************

  function C(n,k)
    integer,intent(in) :: n,k
    real               :: C
    integer :: i

    if (k.eq.0) then
       C=1.0
    else
       C=real(n)/real(k)

       do i=1,k-1
          C=C*real(n-i)/real(i)
       end do
    end if
  end function C

  !****************FONCTIONS BERNSTEIN********************************************

  function b_basis(i,j,x)
    integer,intent(in) :: i,j
    real,   intent(in) :: x
    real               :: b_basis

    b_basis=C(i+j,i)*x**i*(1.0-x)**j

    if ((x.gt.0.0).and.(b_basis.gt.1.0)) then
       print*,'error',b_basis,i,j,x,(1.0-x),C(i+j,i),x**j,(1.0-x)**i
    end if
  end function b_basis

  subroutine init_basis_b(ordre)
    integer,intent(in) :: ordre
    integer :: i
    real,dimension(ordre+1) :: coef

    allocate(base_b(ordre+1,ordre+1))

    coef=0.0
    coef(1)=1.0
    base_b(1,:)=coef

    do i=2,ordre+1
       coef(i-1)=0.0
       coef(i)=1.0
       base_b(i,:)=coef
    end do
  end subroutine init_basis_b


  subroutine free_basis_b
    deallocate(base_b)
  end subroutine free_basis_b


  function eval_polynom_b(pol,x)
    real,dimension(:),intent(in) :: pol
    real             ,intent(in) :: x
    real                         :: eval_polynom_b
    integer                      :: i,j
    real,dimension(:),allocatable:: barycentre1,barycentre2

    eval_polynom_b=0.0
    allocate(barycentre1(size(pol)))
    barycentre1=pol
    do i=1,size(pol)-1
       allocate(barycentre2(size(barycentre1)-1))
       do j=1,size(barycentre2)
          barycentre2(j)=barycentre1(j)*(1-x)+barycentre1(j+1)*(x)
       end do
       deallocate(barycentre1)
       allocate(barycentre1(size(barycentre2)))
       barycentre1=barycentre2
       deallocate(barycentre2)
    end do
    eval_polynom_b=barycentre1(1)

    ! do i=1,size(pol)
    !    eval_polynom_b=eval_polynom_b   &
    !         +pol(i)*b_basis(i-1,size(pol)-i,x)
    ! end do     
  end function eval_polynom_b

  !****************FONCTIONS LAGRANGE********************************************

  function l_basis(ordre,n,x)
    integer,intent(in)       :: ordre
    integer,intent(in)       :: n
    real   ,intent(in)       :: x
    real                     :: l_basis
    real,dimension(ordre+1)  :: xi
    integer                  :: i

    do i=1,ordre+1
       xi(i)=(i-1)*(1.0/ordre)
    end do

    l_basis=1.0

    do i=1,n-1
       l_basis=l_basis*(x-xi(i))/(xi(n)-xi(i))
    end do

    do i=n+1,ordre+1
       l_basis=l_basis*(x-xi(i))/(xi(n)-xi(i))
    end do
  end function l_basis


  subroutine init_basis_l(ordre)
    integer,intent(in) :: ordre
    integer :: i
    real,dimension(ordre+1) :: coef

    allocate(base_l(ordre+1,ordre+1))
    coef=0.0
    coef(1)=1.0
    base_l(1,:)=coef

    do i=2,ordre+1
       coef(i-1)=0.0
       coef(i)=1.0
       base_l(i,:)=coef
    end do
  end subroutine init_basis_l


  subroutine free_basis_l
    deallocate(base_l)
  end subroutine free_basis_l


  function eval_polynom_l(pol,x)
    real,dimension(:),intent(in) :: pol
    real             ,intent(in) :: x
    real                         :: eval_polynom_l
    integer                      :: i

    eval_polynom_l=0.0    
    do i=1,size(pol)
       eval_polynom_l=eval_polynom_l+pol(i)*l_basis(size(pol)-1,i,x)
    end do
  end function eval_polynom_l


  !****************FONCTIONS PASSAGES*******************************************

  subroutine create_B2L
    integer :: i,j

    allocate(B2L(size(base_b,1),size(base_b,1)))

    do j=1,size(base_b,1)
       do i=1,size(base_b,1)
          ! B2L(j,i)=real((j-1))/(size(base_b,1)-1)
          B2L(j,i)=eval_polynom_b(base_b(i,:),real((j-1))/(size(base_b,1)-1))
       end do
    end do
  end subroutine create_B2L


  subroutine free_B2L
    deallocate(B2L)
  end subroutine free_B2L


  subroutine Bernstein2Lagrange(bpol,lpol,nb_elem,DoF)
    real,dimension(nb_elem*DoF),intent(in) :: bpol
    real,dimension(nb_elem*DoF),intent(out) :: lpol
    integer                    ,intent(in)  :: nb_elem
    integer                    ,intent(in)  :: DoF
    real,dimension(DoF):: blocpol

    integer :: i,j
    real    :: x
    do i=1,nb_elem
       lpol((i-1)*DoF+1:i*DoF)=matmul(B2L,bpol((i-1)*DoF+1:i*DoF))
    end do
  end subroutine Bernstein2Lagrange


  subroutine create_L2B
    ! L2B=inv(B2L)
    L2B=LU_inv(B2L)
  end subroutine create_L2B


  subroutine free_L2B
    deallocate(L2B)
  end subroutine free_L2B


  subroutine Lagrange2Bernstein(lpol,bpol)
    real,dimension(:),intent(in)  :: lpol
    real,dimension(:),intent(out) :: bpol
    bpol=matmul(L2B,lpol)
  end subroutine Lagrange2Bernstein



  !!****************FONCTIONS DERIVES*********************************************

  subroutine create_derive(ordre)
    integer,intent(in) :: ordre
    integer :: i

    allocate(D0(ordre+1,ordre+1))
    allocate(D1(ordre+1,ordre+1))  

    D0=0.0
    D1=0.0

    do i=2,ordre+1
       D0(i,i)=-real(i-1)
       D0(i-1,i)=-real(ordre+2-i)
    end do

    do i=1,ordre
       D1(i,i)=-real(ordre+1-i)
       D1(i+1,i)=-real(i)
    end do

    call Full2Sparse(D0,D0_sparse)
    call Full2Sparse(D1,D1_sparse)
  end subroutine create_derive

  subroutine free_derive
    deallocate(D0)
    deallocate(D1)
    call free_sparse_matrix(D0_sparse)
    call free_sparse_matrix(D1_sparse)
  end subroutine free_derive

  subroutine deriv_pol_b(bpol,dbpol)
    real,dimension(:),intent(in)  :: bpol
    real,dimension(:),intent(out) :: dbpol
    dbpol=sparse_matmul(D1_sparse,bpol)-sparse_matmul(D0_sparse,bpol)
  end subroutine deriv_pol_b

end module m_polynom
