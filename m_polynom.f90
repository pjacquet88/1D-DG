module m_polynom
  use m_matrix
  implicit none

  real,dimension(5,5),private :: B2L,L2B,D0,D1
  type(sparse_matrix),private :: D0_sparse,D1_sparse


  type t_polynom_b
     real,dimension(:),allocatable :: coef
     integer                       :: ordre
  end type t_polynom_b

  type t_polynom_l
     real,dimension(:),allocatable :: coef
     integer                       :: ordre
  end type t_polynom_l

  type(t_polynom_b),dimension(:),allocatable :: base_b
  type(t_polynom_l),dimension(:),allocatable :: base_l


contains

  function fact(n)
    integer,intent(in) :: n
    integer            :: fact
    integer :: i

    if (n.eq.0) then
       fact=1
    else

       fact=n

       do i=n-1,2,-1
          fact=fact*i
       end do
    end if
  end function fact

  function C(n,k)
    integer,intent(in) :: n,k
    real               :: C
    integer :: i

    if (k.eq.0) then
       C=1.0
    else
       C=n
       do i=n-1,n-k+1,-1
          C=C*i
       end do

       C=C/fact(k)
    end if
  end function C

  subroutine init_basis_b(ordre)
    integer,intent(in) :: ordre
    integer :: i
    real,dimension(ordre) :: coef

    allocate(base_b(ordre))

    coef=0.0
    coef(1)=1.0
    call init_polynom_b(base_b(1),coef)

    do i=2,ordre
       coef(i-1)=0.0
       coef(i)=1.0
       call init_polynom_b(base_b(i),coef)
    end do

  end subroutine init_basis_b

  subroutine init_basis_l(ordre)
    integer,intent(in) :: ordre
    integer :: i
    real,dimension(ordre) :: coef

    allocate(base_l(ordre))

    coef=0.0
    coef(1)=1.0
    call init_polynom_l(base_l(1),coef)

    do i=2,ordre
       coef(i-1)=0.0
       coef(i)=1.0
       call init_polynom_l(base_l(i),coef)
    end do

  end subroutine init_basis_l
  

  function b_basis(i,j,x)
    integer,intent(in) :: i,j
    real,   intent(in) :: x
    real               :: b_basis

    b_basis=C(i+j,i)*x**j*(1-x)**i

  end function b_basis

  subroutine init_polynom_b(pol,coef)
    type(t_polynom_b),intent(inout) :: pol
    real,dimension(:),intent(in)    :: coef

    pol%ordre=size(coef)-1
    allocate(pol%coef(pol%ordre))
    pol%coef=coef
  end subroutine init_polynom_b

  subroutine free_polynom_b(pol)
    type(t_polynom_b) :: pol

    deallocate(pol%coef)
  end subroutine free_polynom_b

  subroutine free_polynom_l(pol)
    type(t_polynom_l) :: pol

    deallocate(pol%coef)
  end subroutine free_polynom_l
  
  function eval_polynom_b(pol,x)
    type(t_polynom_b),intent(in) :: pol
    real             ,intent(in) :: x
    real                         :: eval_polynom_b

    eval_polynom_b=pol%coef(1)*b_basis(0,4,x)    &
             +pol%coef(2)*b_basis(1,3,x)   &
             +pol%coef(3)*b_basis(2,2,x)   &
             +pol%coef(4)*b_basis(3,1,x)   &
             +pol%coef(5)*b_basis(4,0,x)
  end function eval_polynom_b



  subroutine init_polynom_l(pol,coef)
    type(t_polynom_l),intent(inout) :: pol
    real,dimension(:),intent(in)    :: coef

    pol%ordre=size(coef)-1
    allocate(pol%coef(pol%ordre))
    pol%coef=coef
  end subroutine init_polynom_l

  function l_basis(ordre,n,x)
    integer,intent(in) :: ordre
    integer,intent(in) :: n
    real   ,intent(in) :: x
    real               :: l_basis
    real,dimension(5)  :: xi
    integer            :: i

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

  function eval_polynom_l(pol,x)
    type(t_polynom_l),intent(in) :: pol
    real             ,intent(in) :: x
    real                         :: eval_polynom_l
    integer                      :: i

    eval_polynom_l=0.0
    
    do i=1,pol%ordre+1
       eval_polynom_l=eval_polynom_l+pol%coef(i)*l_basis(pol%ordre,i,x)
    end do
  end function eval_polynom_l

  subroutine create_B2L
    integer :: j
    type(t_polynom_b)::bpol1,bpol2,bpol3,bpol4,bpol5
   
  call init_polynom_b(bpol1,(/1.0,0.0,0.0,0.0,0.0/))
  call init_polynom_b(bpol2,(/0.0,1.0,0.0,0.0,0.0/))
  call init_polynom_b(bpol3,(/0.0,0.0,1.0,0.0,0.0/))
  call init_polynom_b(bpol4,(/0.0,0.0,0.0,1.0,0.0/))
  call init_polynom_b(bpol5,(/0.0,0.0,0.0,0.0,1.0/))

    do j=1,5
     B2L(j,1)=eval_polynom_b(bpol1,(j-1)*0.25)
     B2L(j,2)=eval_polynom_b(bpol2,(j-1)*0.25)
     B2L(j,3)=eval_polynom_b(bpol3,(j-1)*0.25)
     B2L(j,4)=eval_polynom_b(bpol4,(j-1)*0.25)
     B2L(j,5)=eval_polynom_b(bpol5,(j-1)*0.25)
  end do

end subroutine create_B2L

subroutine Bernstein2Lagrange(bpol,lpol)
  type(t_polynom_b),intent(in) :: bpol
  type(t_polynom_l),intent(out) :: lpol

  lpol%coef=matmul(B2L,bpol%coef)
end subroutine Bernstein2Lagrange


subroutine create_L2B
  L2B=inv5(B2L)
end subroutine create_L2B

subroutine create_derive
  integer :: i
  real,dimension(5,5) :: test
  
  D0=0.0
  D1=0.0
  do i=2,5
     D0(i,i)=real(i-1)
     D0(i-1,i)=real(6-i)
  end do

  do i=1,4
     D1(i,i)=real(5-i)
     D1(i+1,i)=real(i)
  end do

  call Full2Sparse(D0,D0_sparse)
  call Full2Sparse(D1,D1_sparse)
  print*,'D0_sparse'
  call print_sparse_matrix(D0_sparse)
  print*,'D1_sparse'
  call print_sparse_matrix(D1_sparse)

end subroutine create_derive

subroutine Lagrange2Bernstein(lpol,bpol)
  type(t_polynom_l),intent(in)  :: lpol
  type(t_polynom_b),intent(out) :: bpol

  bpol%coef=matmul(L2B,lpol%coef)

end subroutine Lagrange2Bernstein

subroutine deriv_pol_b(bpol,dbpol)
  type(t_polynom_b),intent(in)  :: bpol
  type(t_polynom_b),intent(out) :: dbpol
  real,dimension(bpol%ordre+1)  :: coef
  
  coef=sparse_matmul(D0_sparse,bpol%coef)!-sparse_matmul(D0_sparse,bpol%coef)
  call init_polynom_b(dbpol,coef)

end subroutine deriv_pol_b


end module m_polynom
