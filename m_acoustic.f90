module m_acoustic
  use m_polynom
  use m_matrix

  implicit none

  type element
     real :: length
     real :: x_ini
     real :: x_end
     real,dimension(:),allocatable :: coef
     integer :: DoF
     logical :: bernstein
  end type element
  
  real,dimension(:,:),allocatable :: m_loc
  real,dimension(:,:),allocatable :: m_loc_l,m_loc_b
  real,dimension(:,:),allocatable :: m_inv_loc_l,m_inv_loc_b

  real,dimension(15)              :: xx,weight

  real,parameter :: PI=acos(-1.0)

  public  :: m_loc,                                                             &
             init_quadrature,                                                   &
             init_m_loc_b,init_m_loc_l,                                         &
             signal_ini

  private :: xx,weight,                                                         &
             Int_bxb,Int_lxl,                                                   &
             PI


contains



  !****************** QUADRATURE ***********************************************
  subroutine init_quadrature
    integer :: i
    xx(1)=6.0037409897572858D-3
    xx(2)=3.1363303799647048D-2
    xx(3)=7.5896708294786392D-2
    xx(4)=1.3779113431991498D-1
    xx(5)=2.1451391369573058D-1
    xx(6)=3.0292432646121832D-1
    xx(7)=3.9940295300128274D-1
    xx(8)=5.D-1

    weight(1)=1.5376620998058634D-2
    weight(2)=3.5183023744054062D-2
    weight(3)=5.3579610233585968D-2
    weight(4)=6.9785338963077157D-2
    weight(5)=8.3134602908496967D-2
    weight(6)=9.3080500007781106D-2
    weight(7)=9.9215742663555788D-2
    weight(8)=1.0128912096278064D-1

    do i=9,15
       xx(i)=1.0-xx(16-i)
       weight(i)=weight(16-i)
    end do
  end subroutine init_quadrature


  function Int_bxb(bpol1,bpol2)
    type(t_polynom_b),intent(in) :: bpol1
    type(t_polynom_b),intent(in) :: bpol2
    real                         :: Int_bxb
    integer                      :: i

    Int_bxb=0.0

    do i=1,15
       Int_bxb=Int_bxb                  &
            +eval_polynom_b(bpol1,xx(i))&
            *eval_polynom_b(bpol2,xx(i))&
            *weight(i)
    end do
  end function Int_bxb

  function Int_lxl(lpol1,lpol2)
    type(t_polynom_l),intent(in) :: lpol1
    type(t_polynom_l),intent(in) :: lpol2
    real                         :: Int_lxl
    integer                      :: i

    Int_lxl=0.0

    do i=1,15
       Int_lxl=Int_lxl                  &
            +eval_polynom_l(lpol1,xx(i))&
            *eval_polynom_l(lpol2,xx(i))&
            *weight(i)
    end do
  end function Int_lxl



  !**************** MASS MATRIX ************************************************
  subroutine init_m_loc_b(DoF)
    integer,intent(in) :: DoF
    integer            :: i,j

    allocate(m_loc(DoF,DoF))

    ! do i=1,DoF
    !     m_loc(i,i)=Int_bxb(base_b(i),base_b(i))
    !    do j=i+1,DoF
    !       m_loc(i,j)=Int_bxb(base_b(i),base_b(j))
    !       m_loc(i,DoF-j+1)=m_loc(i,j)
    !    end do
    ! end do

    do i=1,DoF
       do j=1,DoF
          m_loc(i,j)=Int_bxb(base_b(i),base_b(j))
       end do
    end do

  end subroutine init_m_loc_b

  subroutine init_m_loc_l(DoF)
    integer,intent(in) :: DoF
    integer            :: i,j

    allocate(m_loc(DoF,DoF))

    ! do i=1,DoF
    !     m_loc(i,i)=Int_lxl(base_l(i),base_l(i))
    !    do j=i+1,DoF
    !       m_loc(i,j)=Int_lxl(base_l(i),base_l(j))
    !       m_loc(i,DoF-j+1)=m_loc(i,j)
    !    end do
    ! end do

    do i=1,DoF
       do j=1,DoF
          m_loc(i,j)=Int_lxl(base_l(i),base_l(j))
       end do
    end do
  end subroutine init_m_loc_l

  !******************* INIT PROBLEM **************************************

  function signal_ini(x,a)
    real,intent(in) :: x
    character(len=*),intent(in) :: a
    real                        :: signal_ini
    
    real                        :: f

    f=25.0

    if (a.eq.'sinus') then
       signal_ini=sin(4*PI*x)
    else
       signal_ini=(x-0.5)*exp(-(2.0*PI*(x-0.5)*2.0)**2.0)
    end if
  end function signal_ini


  subroutine init_element(elem,nb_elem,DoF,signal,total_length,bernstein)
    type(element),dimension(nb_elem),intent(inout) :: elem
    integer                         ,intent(in)    :: nb_elem
    integer                         ,intent(in)    :: DoF
    logical                         ,intent(in)    :: bernstein
    character(len=*)                ,intent(in)    :: signal
    real                            ,intent(in)    :: total_length

    integer :: i,j
    real    :: x

    elem(1)%length=total_length/nb_elem
    elem(1)%x_ini=0.0
    elem(1)%x_end=elem(1)%length
    elem(1)%DoF=DoF
    allocate(elem(1)%coef(DoF))
    elem(1)%bernstein=bernstein
    
    do j=1,elem(1)%DoF
       x=elem(1)%x_ini+(j-1)*elem(1)%length/(elem(1)%DoF-1)
       elem(1)%coef(j)=signal_ini(x,signal)
    end do

    do i=2,nb_elem
       elem(i)%length=elem(i-1)%length
       elem(i)%x_ini=elem(i-1)%x_end
       elem(i)%x_end=elem(i)%x_ini+elem(i)%length
       elem(i)%DoF=elem(i-1)%DoF
       allocate(elem(i)%coef(DoF))
       elem(i)%bernstein=bernstein
       
       do j=1,elem(i)%DoF
          x=elem(i)%x_ini+(j-1)*elem(i)%length/(elem(i)%DoF-1)
          elem(i)%coef(j)=signal_ini(x,signal)
       end do
    end do

    print*,'test',signal
    
  end subroutine init_element

  subroutine print_sol(elem,N)
    type(element),dimension(:),intent(in) :: elem
    integer                   ,intent(in) :: N
    real,dimension(:),allocatable         :: coef

    integer :: i,j
    character(len=20) :: F_NAME

    
    write(F_NAME,"(A,I0,'.dat')") "sol",N

    open(unit=2, file=F_NAME, action="write")


    
    do i=1,size(elem)
       if (elem(i)%bernstein) then
          allocate(coef(elem(i)%DoF))
          coef=matmul(L2B,elem(i)%coef)
          do j=1,elem(i)%DoF
             write(2,*),elem(i)%x_ini+(j-1)*elem(i)%length/(elem(i)%DoF-1),coef(j)
          end do
          deallocate(coef)
       else
          do j=1,elem(i)%DoF
             write(2,*),elem(i)%x_ini+(j-1)*elem(i)%length/(elem(i)%DoF-1),elem(i)%coef(j)
          end do
       end if
    end do
  end subroutine print_sol
  
end module m_acoustic
