! This module contains matrices and functions used in Lagrange and Bernstein
! polynomial definitions
module m_polynom
  use m_kind
  use m_matrix
  implicit none

  real(mp),dimension(:,:),allocatable :: B2L,L2B,D0,D1
  real(mp),dimension(:,:),allocatable :: base_b
  real(mp),dimension(:,:),allocatable :: base_l
  type(sparse_matrix)                 :: D0_sparse,D1_sparse


  public  :: eval_polynom_b,eval_polynom_l,                                     &
             init_basis_b,init_basis_l,                                         &
             create_L2B,create_B2L,create_derive,                               &
             free_L2B,free_B2L,free_derive,                                     &
             Bernstein2Lagrange,Lagrange2Bernstein,B2L,L2B,D1,D0

  private :: D0_sparse,D1_sparse,                                               &
             C,b_basis,l_basis
contains


  ! Calculates the combination C(n,k)
  function C(n,k)
    integer,intent(in) :: n,k
    real(mp)           :: C
    integer            :: i

    if (k.eq.0) then
       C=1.0_mp
    else
       C=real(n)/real(k)

       do i=1,k-1
          C=C*real(n-i)/real(i)
       end do
    end if
  end function C


  ! Evaluates one Bernstein polynomial basis B(i,j) at x, with x in [0,1]
  function b_basis(i,j,x)
    integer ,intent(in) :: i,j
    real(mp),intent(in) :: x
    real(mp)            :: b_basis

    b_basis=C(i+j,i)*x**i*(1.0_mp-x)**j

    if ((x.gt.0.0_mp).and.(b_basis.gt.1.0_mp)) then
       print*,'error',b_basis,i,j,x,(1.0_mp-x),C(i+j,i),x**j,(1.0_mp-x)**i
    end if
  end function b_basis

  ! Initializes the Bernstein basis
  subroutine init_basis_b(ordre)
    integer,intent(in) :: ordre

    integer                     :: i
    real(mp),dimension(ordre+1) :: coef

    allocate(base_b(ordre+1,ordre+1))

    coef=0.0_mp
    coef(1)=1.0_mp
    base_b(1,:)=coef

    do i=2,ordre+1
       coef(i-1)=0.0_mp
       coef(i)=1.0_mp
       base_b(i,:)=coef
    end do
  end subroutine init_basis_b

  ! Deallocates the Bernstein basis
  subroutine free_basis_b
    deallocate(base_b)
  end subroutine free_basis_b


  ! Evaluates a Bernstein polynomial using De Calsteljau algorithm
  function eval_polynom_b(pol,x)
    real(mp),dimension(:),intent(in)  :: pol
    real(mp)             ,intent(in)  :: x
    real(mp)                          :: eval_polynom_b
    real(mp),dimension(:),allocatable :: barycentre1,barycentre2
    integer                           :: i,j

    eval_polynom_b=0.0_mp
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
  end function eval_polynom_b


  ! Evaluates one Lagrange polynomial basis at x, with x in [O,1]
  function l_basis(ordre,n,x)
    integer ,intent(in)       :: ordre
    integer ,intent(in)       :: n
    real(mp),intent(in)       :: x

    real(mp)                    :: l_basis
    real(mp),dimension(ordre+1) :: xi
    integer                     :: i

    do i=1,ordre+1
       xi(i)=(i-1)*(1.0_mp/ordre)
    end do

    l_basis=1.0_mp

    do i=1,n-1
       l_basis=l_basis*(x-xi(i))/(xi(n)-xi(i))
    end do

    do i=n+1,ordre+1
       l_basis=l_basis*(x-xi(i))/(xi(n)-xi(i))
    end do
  end function l_basis


  ! Initializes the Lagrange polynomial basis
  subroutine init_basis_l(ordre)
    integer,intent(in) :: ordre
    integer :: i
    real(mp),dimension(ordre+1) :: coef

    allocate(base_l(ordre+1,ordre+1))
    coef=0.0_mp
    coef(1)=1.0_mp
    base_l(1,:)=coef

    do i=2,ordre+1
       coef(i-1)=0.0_mp
       coef(i)=1.0_mp
       base_l(i,:)=coef
    end do
  end subroutine init_basis_l


  ! Deallocates the Lagrange basis
  subroutine free_basis_l
    deallocate(base_l)
  end subroutine free_basis_l


  ! Evaluates a Lagrange polynomial at x, with x in [0,1]
  function eval_polynom_l(pol,x)
    real(mp),dimension(:),intent(in) :: pol
    real(mp)             ,intent(in) :: x
    real(mp)                         :: eval_polynom_l
    integer                          :: i

    eval_polynom_l=0.0_mp
    do i=1,size(pol)
       eval_polynom_l=eval_polynom_l+pol(i)*l_basis(size(pol)-1,i,x)
    end do
  end function eval_polynom_l


  !**************** LAGRANGE <--> BERNSTEIN *************************************

  ! Creates a matrix that transforms a vector of Bernstein polynomial coefficients
  ! into a vector of Lagrange polynomial coefficients
  subroutine create_B2L
    integer  :: i,j
    real(mp) :: x

    allocate(B2L(size(base_b,1),size(base_b,1)))

    do j=1,size(base_b,1)
       x=real((j-1))/(size(base_b,1)-1)
       do i=1,size(base_b,1)
          B2L(j,i)=eval_polynom_b(base_b(i,:),x)
       end do
    end do
  end subroutine create_B2L


  ! Deallocates the B2L matrix
  subroutine free_B2L
    deallocate(B2L)
  end subroutine free_B2L


  ! Subroutine which changes a Bernstein polynomial coefficients into
  ! a Lagrangian polynomial coefficients
  subroutine Bernstein2Lagrange(bpol,lpol,nb_elem,DoF)
    real(mp),dimension(nb_elem*DoF),intent(in)  :: bpol
    real(mp),dimension(nb_elem*DoF),intent(out) :: lpol
    integer                        ,intent(in)  :: nb_elem
    integer                        ,intent(in)  :: DoF
    real(mp),dimension(DoF):: blocpol

    integer  :: i,j
    real(mp) :: x

    do i=1,nb_elem
       lpol((i-1)*DoF+1:i*DoF)=matmul(B2L,bpol((i-1)*DoF+1:i*DoF))
    end do
  end subroutine Bernstein2Lagrange


  ! Creates a matrix that transforms a vector of Lagrange polynomial coefficients
  ! into a vector of Bernstein polynomial coefficients
  subroutine create_L2B
    L2B=LU_inv(B2L)
  end subroutine create_L2B


  ! Deallocates the L2B matrix
  subroutine free_L2B
    deallocate(L2B)
  end subroutine free_L2B



  ! Subroutine which changes a Lagrange polynomial coefficients into
  ! a Bernstein polynomial coefficients
  subroutine Lagrange2Bernstein(lpol,bpol)
    real(mp),dimension(:),intent(in)  :: lpol
    real(mp),dimension(:),intent(out) :: bpol
    bpol=matmul(L2B,lpol)
  end subroutine Lagrange2Bernstein


  !!****************FONCTIONS DERIVES*******************************************

  ! Creates all the Bernstein derivative matrix
  ! D0 = derivative matrix by lambda_0
  ! D1 = derivative matrix by lambda_1
  subroutine create_derive(ordre)
    integer,intent(in) :: ordre
    integer :: i

    allocate(D0(ordre+1,ordre+1))
    allocate(D1(ordre+1,ordre+1))  

    D0=0.0_mp
    D1=0.0_mp

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


  ! Deallocates the Bernstein derivative matrix
  subroutine free_derive
    deallocate(D0)
    deallocate(D1)
    call free_sparse_matrix(D0_sparse)
    call free_sparse_matrix(D1_sparse)
  end subroutine free_derive

  ! Gets the global sparse derivative matrix
  subroutine deriv_pol_b(bpol,dbpol)
    real(mp),dimension(:),intent(in)  :: bpol
    real(mp),dimension(:),intent(out) :: dbpol
    dbpol=sparse_matmul(D1_sparse,bpol)-sparse_matmul(D0_sparse,bpol)
  end subroutine deriv_pol_b

end module m_polynom
