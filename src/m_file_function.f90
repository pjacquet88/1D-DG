! This module contain some function concerning files and vectors
module m_file_function
  use m_kind
  use m_polynom
  implicit none

  public  :: load,print_vect,laplace_filter
  private :: norm,nb_line

contains

  !*************** PRIVATE ******************************************************
  ! Calculate the norm of a vector U
  function norm(U)
    real(mp),dimension(:),intent(in) :: U
    real(mp)                         :: norm

    norm=sqrt(dot_product(U,U))
  end function norm



  ! Get the number of line of a file
  function nb_line(fichier)
    character(len=*),intent(in) :: fichier
    integer :: ierr ! successful (if == 0) ou failed (si /= 0)
    integer :: nb_line
    open(12,file=fichier, form="formatted", iostat=ierr,status="old")

    nb_line = 0
    do while (ierr==0)
       read(12,*,iostat=ierr)
       if (ierr==0) nb_line = nb_line + 1
    enddo
    close(12)
  end function nb_line



  !*************** PUBLIC   *****************************************************
  ! Load a vector in a specific file .dat
  subroutine load(U,char,N)
    real(mp),dimension(:),intent(inout):: U
    character(len=*) ,intent(in)       :: char
    integer          ,intent(in)       :: N
    real(mp)                           :: dummy
    integer                            :: i
    character(len=50)                  :: fichier

    U=0.0_mp
    write(fichier,"(A,A,I0,'.dat')") "../Bernstein/Files/",char,N
    
    open(unit=456,file=fichier, form="formatted",action='read')
    do i=1,size(U)
       read(456,*) dummy,U(i)
    end do
    close(456)
  end subroutine load


   !***************** SORTIE DE RESULTAT *********************************

  ! Print a vector in a specific file
  subroutine print_vect(vector,nb_elem,DoF,dx,bernstein,N,name)
    real(mp),dimension(nb_elem*DoF),intent(in) :: vector
    integer                        ,intent(in) :: nb_elem
    integer                        ,intent(in) :: DoF
    real(mp)                       ,intent(in) :: dx
    logical                        ,intent(in) :: bernstein
    integer                        ,intent(in) :: N
    character(len=*)               ,intent(in) :: name
    
    integer :: i,j
    character(len=20)           :: F_NAME
    real(mp)                    :: x,ddx
    real(mp),dimension(DoF)     :: vector_work

    ddx=dx/(DoF-1)
    x=0.0_mp
    
    if (N.eq.1300) then
       print*,'test kind vector', kind(vector)
       print*,'vector(1) :',vector(1)
       STOP
    end if
    
    write(F_NAME,"(A,A,I0,'.dat')") "../Files/",name,N
    open(unit=2, file=F_NAME, action="write")
    if (bernstein) then
       do i=1,nb_elem-1
          vector_work=matmul(B2L,vector(Dof*(i-1)+1:Dof*i))
          do j=1,DoF-1
             x=(i-1)*dx+(j-1)*ddx
             write(2,*)x,vector_work(j)
          end do
       end do
       vector_work=matmul(B2L,vector(Dof*(nb_elem-1)+1:Dof*nb_elem))
       do j=1,DoF
          x=(nb_elem-1)*dx+(j-1)*ddx
          write(2,*)x,vector_work(j)
       end do
    else
       do i=1,nb_elem-1
          do j=1,DoF-1
             x=(i-1)*dx+(j-1)*ddx
             write(2,*)x,vector((i-1)*DoF+j)
          end do
       end do
       do j=1,DoF
          x=(nb_elem-1)*dx+(j-1)*ddx
          write(2,*)x,vector((nb_elem-1)*DoF+j)
       end do

    end if
    close(2)
  end subroutine print_vect


  ! Add a laplacian filter on a vector
  function laplace_filter(u)
    real(mp),dimension(:),intent(in)    :: u
    real(mp),dimension(size(u),size(u)) :: L,Linv,test
    integer                             :: i,j
    real(mp),dimension(size(u))         :: laplace_filter

    L=0.0_mp

    do i=30,size(u)-30
       L(i,i)=2
       L(i,i-1)=-1
       L(i,i+1)=-1
    end do

    L(size(u),size(u))=2
    L(size(u),size(u)-1)=-1

    laplace_filter=matmul(L,u)
  end function laplace_filter

end module m_file_function
