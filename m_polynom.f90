module m_polynom
  use m_matrix
  implicit none

  real               ,dimension(:,:),allocatable :: B2L,L2B,D0,D1
  type(sparse_matrix)                            :: D0_sparse,D1_sparse

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
  

  public  :: t_polynom_b,t_polynom_l,base_b,base_l,                             &
             init_polynom_b,init_polynom_l,                                     &
             free_polynom_b,free_polynom_l,                                     &
             eval_polynom_b,eval_polynom_l,                                     &
             init_basis_b,init_basis_l,                                         &
             create_L2B,create_B2L,create_derive,                               &
             Bernstein2Lagrange,Lagrange2Bernstein,B2L
  
  private :: L2B,D0,D1,D0_sparse,D1_sparse,                                 &
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
  
  subroutine init_polynom_b(pol,coef)
    type(t_polynom_b),intent(inout) :: pol
    real,dimension(:),intent(in)    :: coef

    pol%ordre=size(coef)-1
    allocate(pol%coef(pol%ordre))
    pol%coef=coef
  end subroutine init_polynom_b
   
  function b_basis(i,j,x)
    integer,intent(in) :: i,j
    real,   intent(in) :: x
    real               :: b_basis

    b_basis=C(i+j,i)*x**j*(1.0-x)**i

    if ((x.gt.0.0).and.(b_basis.gt.1.0)) then
       print*,'error',b_basis,i,j,x,(1.0-x),C(i+j,i),x**j,(1.0-x)**i
    end if
  end function b_basis

   subroutine init_basis_b(ordre)
    integer,intent(in) :: ordre
    integer :: i
    real,dimension(ordre+1) :: coef

    allocate(base_b(ordre+1))

    coef=0.0
    coef(1)=1.0
    call init_polynom_b(base_b(1),coef)

    do i=2,ordre+1
       coef(i-1)=0.0
       coef(i)=1.0
       call init_polynom_b(base_b(i),coef)
    end do

  end subroutine init_basis_b  
  
  function eval_polynom_b(pol,x)
    type(t_polynom_b),intent(in) :: pol
    real             ,intent(in) :: x
    real                         :: eval_polynom_b
    integer                      :: i

    eval_polynom_b=0.0
    
    do i=1,pol%ordre+1
       eval_polynom_b=eval_polynom_b   &
            +pol%coef(i)*b_basis(i-1,pol%ordre-i+1,x)
    end do
       
  end function eval_polynom_b

  
  subroutine free_polynom_b(pol)
    type(t_polynom_b) :: pol

    deallocate(pol%coef)
  end subroutine free_polynom_b



  
  !****************FONCTIONS LAGRANGE********************************************
  
  subroutine init_polynom_l(pol,coef)
    type(t_polynom_l),intent(inout) :: pol
    real,dimension(:),intent(in)    :: coef

    pol%ordre=size(coef)-1
    allocate(pol%coef(pol%ordre))
    pol%coef=coef
  end subroutine init_polynom_l

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

    allocate(base_l(ordre+1))
    coef=0.0
    coef(1)=1.0
    call init_polynom_l(base_l(1),coef)

    do i=2,ordre+1
       coef(i-1)=0.0
       coef(i)=1.0
       call init_polynom_l(base_l(i),coef)
    end do

  end subroutine init_basis_l

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

  subroutine free_polynom_l(pol)
    type(t_polynom_l) :: pol

    deallocate(pol%coef)
  end subroutine free_polynom_l



  
  !****************FONCTIONS PASSAGES*******************************************
  
  subroutine create_B2L
    integer :: i,j

    allocate(B2L(size(base_b),size(base_b)))

    do j=1,size(base_b)
       do i=1,size(base_b)
          B2L(j,i)=eval_polynom_b(base_b(i),real((j-1))/(size(base_b)-1))
       end do
    end do
    
end subroutine create_B2L

subroutine Bernstein2Lagrange(bpol,lpol)
  type(t_polynom_b),intent(in) :: bpol
  type(t_polynom_l),intent(out) :: lpol

  call init_polynom_l(lpol,matmul(B2L,bpol%coef))

end subroutine Bernstein2Lagrange


subroutine create_L2B
! L2B=inv(B2L)
  L2B=LU_inv(B2L)
end subroutine create_L2B

subroutine Lagrange2Bernstein(lpol,bpol)
  type(t_polynom_l),intent(in)  :: lpol
  type(t_polynom_b),intent(out) :: bpol

  call init_polynom_b(bpol,matmul(L2B,lpol%coef))

end subroutine Lagrange2Bernstein



!!****************FONCTIONS DERIVES*********************************************

subroutine create_derive(ordre)
  integer,intent(in) :: ordre
  integer :: i
  real,dimension(ordre+1,ordre+1) :: test

  allocate(D0(ordre+1,ordre+1))
  allocate(D1(ordre+1,ordre+1))  
  
  D0=0.0
  D1=0.0
  
  do i=2,ordre+1
     D0(i,i)=real(i-1)
     D0(i-1,i)=real(ordre+2-i)
  end do

  do i=1,ordre
     D1(i,i)=real(ordre+1-i)
     D1(i+1,i)=real(i)
  end do

  test=D1-D0

  ! print*,'test'
  ! do i=1,ordre+1
  !    print*,test(i,:)
  ! end do

  call Full2Sparse(D0,D0_sparse)
  call Full2Sparse(D1,D1_sparse)

end subroutine create_derive


subroutine deriv_pol_b(bpol,dbpol)
  type(t_polynom_b),intent(in)  :: bpol
  type(t_polynom_b),intent(out) :: dbpol
  real,dimension(bpol%ordre+1)  :: coef
  
  !coef=matmul(D1,bpol%coef)-matmul(D0,bpol%coef)
  coef=sparse_matmul(D1_sparse,bpol%coef)-sparse_matmul(D0_sparse,bpol%coef)
  call init_polynom_b(dbpol,coef)

end subroutine deriv_pol_b


end module m_polynom
