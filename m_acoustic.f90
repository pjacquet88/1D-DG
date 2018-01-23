module m_acoustic
  use m_polynom
  use m_matrix

  implicit none

  type acoustic_problem
     integer                       :: DoF
     integer                       :: nb_elem
     logical                       :: bernstein
     logical                       :: F_forte
     real                          :: dx
     real                          :: dt
     real                          :: total_length
     real                          :: final_time
     real                          :: alpha
     real,dimension(:),allocatable :: U
     real,dimension(:),allocatable :: P
     character(len=20)             :: boundaries
     character(len=20)             :: signal
     type(sparse_matrix)           :: Ap
     type(sparse_matrix)           :: Av
  end type acoustic_problem

  real,dimension(15)              :: xx,weight

  real,parameter :: PI=acos(-1.0)

  public  :: init_problem,one_time_step,init_ApAv,                              &
             signal_ini,                                                        &
             print_vect,print_sol

  private :: xx,weight,                                                         &
             init_quadrature,Int_bxb,Int_lxl, init_m_loc,init_s_loc,            &
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


    !******************* INIT PROBLEM **************************************
  
  function signal_ini(x,dt,a)
    real,intent(in) :: x
    real,intent(in) :: dt
    character(len=*),intent(in) :: a
    real                        :: signal_ini

    real                        :: xt

    xt=x+dt
    
    if (a.eq.'sinus') then
       signal_ini=sin(4*PI*xt)
       ! signal_ini=sin(xt)
    else if (a.eq.'creneau') then
       if ((x.lt.0.55).and.(x.gt.0.45)) then
          signal_ini=0.5
       else
          signal_ini=0.0
       end if
    else if (a.eq.'x') then
       signal_ini=xt**5
    else
       signal_ini=(xt-0.5)*exp(-(2.0*PI*(xt-0.5)*2.0)**2.0)*20
    end if
    
end function signal_ini


  subroutine init_problem(problem,nb_elem,DoF,total_length,final_time,alpha,    &
                          bernstein,F_forte,signal,boundaries)
    type(acoustic_problem),intent(out) :: problem
    integer               ,intent(in)  :: nb_elem
    integer               ,intent(in)  :: DoF
    real                  ,intent(in)  :: total_length
    real                  ,intent(in)  :: final_time
    real                  ,intent(in)  :: alpha
    logical               ,intent(in)  :: bernstein
    logical               ,intent(in)  :: F_forte
    character(len=*)      ,intent(in)  :: signal
    character(len=*)      ,intent(in)  :: boundaries
    
    integer :: i,j
    real    :: x,ddx

    problem%DoF=DoF
    problem%nb_elem=nb_elem
    problem%total_length=total_length
    problem%dx=total_length/(nb_elem)
    problem%alpha=alpha
    problem%bernstein=bernstein
    problem%F_forte=F_forte
    problem%boundaries=boundaries

    problem%dt=2.05*(0.0682*problem%dx)/10.0

    allocate(problem%U(DoF*nb_elem),problem%P(DoF*nb_elem))
    problem%U=0.0
    problem%P=0.0

    ddx=problem%dx/(DoF-1)
    x=0.0

    do i=1,nb_elem
       do j=1,DoF
          x=(i-1)*problem%dx+(j-1)*ddx
          problem%U((i-1)*DoF+j)=signal_ini(x,0.0,signal)
          problem%P((i-1)*DoF+j)=signal_ini(x,-problem%dt,signal)
       end do
    end do

    if (bernstein) then
       do i=1,nb_elem
          problem%U(Dof*(i-1)+1:Dof*i)=matmul(L2B,problem%U(Dof*(i-1)+1:Dof*i))
          problem%P(Dof*(i-1)+1:Dof*i)=matmul(L2B,problem%P(Dof*(i-1)+1:Dof*i))
       end do
    end if

    print*,'nb_elem              ::',nb_elem
    print*,'DoF                  ::',DoF
    print*,'total_length         ::',total_length
    print*,'final_time           ::',final_time
    print*,'dx                   ::',problem%dx
    print*,'dt                   ::',problem%dt
    print*,'Conditions aux bords ::','   ',boundaries
    print*,'Pénalisation alpha   ::',alpha
    print*,' '
    if (bernstein) then
       print*,'Le problème est résolu avec des éléments de Bernstein'
    else
       print*,'Le problème est résolu avec des éléments de Lagrange'
    end if
    if (F_forte) then
       print*,'Le problème est résolu en formulation forte'
    else
       print*,'Le problème est résolu en formulation faible'
    end if
    print*,'------------------------------------------------------------'
    
  end subroutine init_problem
  

  !**************** RESOLUTION DU PROBLEM *************************************

  subroutine init_m_loc(m_loc,m_inv_loc,bernstein)
    real,dimension(:,:),intent(out) :: m_loc
    real,dimension(:,:),intent(out) :: m_inv_loc
    logical            ,intent(in)  :: bernstein
    integer                         :: DoF
    integer :: i,j


    DoF=size(m_loc,1)
    
    if (bernstein) then
       do i=1,DoF
          m_loc(i,i)=Int_bxb(base_b(i),base_b(i))
          do j=i+1,DoF
             m_loc(i,j)=Int_bxb(base_b(i),base_b(j))
             m_loc(j,i)=m_loc(i,j)
          end do
       end do
    else
       do i=1,DoF
          m_loc(i,i)=Int_lxl(base_l(i),base_l(i))
          do j=i+1,DoF
             m_loc(i,j)=Int_lxl(base_l(i),base_l(j))
             m_loc(j,i)=m_loc(i,j)
          end do
       end do
    end if
    
    m_inv_loc=LU_inv(m_loc)

  end subroutine init_m_loc
  
  subroutine init_s_loc(s_loc,bernstein)
    real,dimension(:,:),intent(out) :: s_loc
    logical            ,intent(in)  :: bernstein
    integer                         :: DoF
    integer :: i,j
    type(t_polynom_b) :: bpol,dbpol
    type(t_polynom_l) :: dlpol

    DoF=size(s_loc,1)

    if (bernstein) then
       s_loc=D1-D0
    else
       do i=1,DoF
          do j=1,DoF
             call Lagrange2Bernstein(base_l(j),bpol)
             call deriv_pol_b(bpol,dbpol)
             call Bernstein2Lagrange(dbpol,dlpol)

             s_loc(i,j)=Int_lxl(base_l(i),dlpol)

          end do
       end do
    end if
  end subroutine init_s_loc


  subroutine init_ApAv(problem)
    type(acoustic_problem),intent(inout) :: problem

    real,dimension(problem%nb_elem*problem%DoF,                                 &
                   problem%nb_elem*problem%DoF) :: s_glob,MINV,F1,F2,Av_full,Ap_full
    real,dimension(problem%DoF,problem%DoF)     :: m_loc
    real,dimension(problem%DoF,problem%DoF)     :: m_inv_loc
    real,dimension(problem%DoF,problem%DoF)     :: s_loc
    real                                        :: alpha
    
    integer                 :: i,j
    integer                 :: DoF,nb_elem

    DoF=problem%DoF           ! For sake of lisibility
    nb_elem=problem%nb_elem

    Ap_full=0.0
    Av_full=0.0
    if (problem%boundaries.eq.'periodique') then


       alpha=problem%alpha

       F1=0.0

       F1(1,1)=0.5+alpha
       F1(1,DoF*nb_elem)=-0.5-alpha

       F1(DoF,DoF)=-0.5+alpha
       F1(DoF,DoF+1)=0.5-alpha

       do i=2,nb_elem-1
          F1(Dof*(i-1)+1,Dof*(i-1)+1)=0.5+alpha
          F1(Dof*(i-1)+1,Dof*(i-1))=-0.5-alpha

          F1(Dof*i,Dof*i)=-0.5+alpha
          F1(Dof*i,Dof*i+1)=0.5-alpha
       end do
       F1(Dof*(nb_elem-1)+1,Dof*(nb_elem-1)+1)=0.5+alpha
       F1(Dof*(nb_elem-1)+1,Dof*(nb_elem-1))=-0.5-alpha

       F1(Dof*nb_elem,Dof*nb_elem)=-0.5+alpha
       F1(Dof*nb_elem,1)=0.5-alpha

       alpha=-problem%alpha

       F2=0.0

       F2(1,1)=0.5+alpha
       F2(1,DoF*nb_elem)=-0.5-alpha

       F2(DoF,DoF)=-0.5+alpha
       F2(DoF,DoF+1)=0.5-alpha

       do i=2,nb_elem-1
          F2(Dof*(i-1)+1,Dof*(i-1)+1)=0.5+alpha
          F2(Dof*(i-1)+1,Dof*(i-1))=-0.5-alpha

          F2(Dof*i,Dof*i)=-0.5+alpha
          F2(Dof*i,Dof*i+1)=0.5-alpha
       end do
       F2(Dof*(nb_elem-1)+1,Dof*(nb_elem-1)+1)=0.5+alpha
       F2(Dof*(nb_elem-1)+1,Dof*(nb_elem-1))=-0.5-alpha

       F2(Dof*nb_elem,Dof*nb_elem)=-0.5+alpha
       F2(Dof*nb_elem,1)=0.5-alpha

    else
       print*,'ELSE',problem%boundaries

       
       alpha=problem%alpha

       F1=0.0

       F1(1,1)=1.0!+2.0*alpha
       F1(1,DoF*nb_elem)=0.0

       F1(DoF,DoF)=-0.5+alpha
       F1(DoF,DoF+1)=0.5-alpha

       do i=2,nb_elem-1
          F1(Dof*(i-1)+1,Dof*(i-1)+1)=0.5+alpha
          F1(Dof*(i-1)+1,Dof*(i-1))=-0.5-alpha

          F1(Dof*i,Dof*i)=-0.5+alpha
          F1(Dof*i,Dof*i+1)=0.5-alpha
       end do
       F1(Dof*(nb_elem-1)+1,Dof*(nb_elem-1)+1)=0.5+alpha
       F1(Dof*(nb_elem-1)+1,Dof*(nb_elem-1))=-0.5-alpha

       F1(Dof*nb_elem,Dof*nb_elem)=-1.0!+2.0*alpha
       F1(Dof*nb_elem,1)=0.0

       alpha=-problem%alpha

       F2=0.0

       F2(1,1)=0.0
       F2(1,DoF*nb_elem)=0.0

       F2(DoF,DoF)=-0.5+alpha
       F2(DoF,DoF+1)=0.5-alpha

       do i=2,nb_elem-1
          F2(Dof*(i-1)+1,Dof*(i-1)+1)=0.5+alpha
          F2(Dof*(i-1)+1,Dof*(i-1))=-0.5-alpha

          F2(Dof*i,Dof*i)=-0.5+alpha
          F2(Dof*i,Dof*i+1)=0.5-alpha
       end do
       F2(Dof*(nb_elem-1)+1,Dof*(nb_elem-1)+1)=0.5+alpha
       F2(Dof*(nb_elem-1)+1,Dof*(nb_elem-1))=-0.5-alpha

       F2(Dof*nb_elem,Dof*nb_elem)=0.0
       F2(Dof*nb_elem,1)=0.0
    end if

    call init_quadrature
    call init_m_loc(m_loc,m_inv_loc,problem%bernstein)
    call init_s_loc(s_loc,problem%bernstein)

    MINV=0.0
    s_glob=0.0

    do i=1,nb_elem
       MINV(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=m_inv_loc
    end do
    
    do i=1,nb_elem
       s_glob(DoF*(i-1)+1:Dof*i,DoF*(i-1)+1:Dof*i)=s_loc
    end do

    if (problem%bernstein) then

       if (problem%boundaries.eq.'neumann') then
          print*,'JE SUIS PASSE PAR LA'
          Ap_full=(s_glob+matmul(MINV,F2))
          Av_full=(s_glob+matmul(MINV,F1))
       else
       Ap_full=(s_glob+matmul(MINV,F1))
       Av_full=(s_glob+matmul(MINV,F2))
    end if
    
       if (.not.problem%F_forte) then
          Ap_full=-transpose(Ap_full)
          Av_full=-transpose(Av_full)
       end if

       call Full2Sparse(Ap_full,problem%Ap)
       call Full2Sparse(Av_full,problem%Av)
       problem%Ap%Values=(problem%dt/problem%dx)*problem%Ap%Values
       problem%Av%Values=(problem%dt/problem%dx)*problem%Av%Values
    else

       if (problem%boundaries.eq.'neumann') then
          print*,'JE SUIS PASSE PAR LA'
          Ap_full=s_glob+F2
          Av_full=s_glob+F1
       else
          Ap_full=s_glob+F1
          Av_full=s_glob+F2   
       end if
       
       if (.not.problem%F_forte) then
          Ap_full=-transpose(Ap_full)
          Av_full=-transpose(Av_full)
       end if

       Ap_full=matmul(MINV,Ap_full)
       Av_full=matmul(MINV,Av_full)

       call Full2Sparse(Ap_full,problem%Ap)
       call Full2Sparse(Av_full,problem%Av)
       problem%Ap%Values=(problem%dt/problem%dx)*problem%Ap%Values
       problem%Av%Values=(problem%dt/problem%dx)*problem%Av%Values
       
    end if
    
    print*,'Nombre de termes non nuls de Ap,Av ::',problem%Ap%NNN
    print*,'Taille des matrices Ap,AV          ::',problem%Ap%nb_ligne,'x', &
         problem%Ap%nb_ligne,'=',problem%Ap%nb_ligne**2
    print*,'Ratio :',real(problem%Ap%NNN)/problem%Ap%nb_ligne**2

  end subroutine init_ApAv
  


  subroutine one_time_step(problem)
    type(acoustic_problem),intent(inout) :: problem


       problem%U=problem%U-sparse_matmul(problem%Av,problem%P)
       problem%P=problem%P-sparse_matmul(problem%Ap,problem%U)
    
  end subroutine one_time_step




 !***************** SORTIE DE RESULTAT *********************************
  
  subroutine print_vect(vector,nb_elem,DoF,dx,bernstein,N,name)
    real,dimension(nb_elem*DoF),intent(in) :: vector

    integer                    ,intent(in) :: nb_elem
    integer                    ,intent(in) :: DoF
    real                       ,intent(in) :: dx
    logical                    ,intent(in) :: bernstein
    integer                    ,intent(in) :: N
    character(len=*)           ,intent(in) :: name
    
    integer :: i,j
    character(len=20)           :: F_NAME
    real                        :: x,ddx
    real,dimension(DoF) :: vector_work

    ddx=dx/(DoF-1)
    x=0.0


    write(F_NAME,"(A,A,I0,'.dat')") "fichier/",name,N

    open(unit=2, file=F_NAME, action="write")

    if (bernstein) then
       do i=1,nb_elem
          vector_work=matmul(B2L,vector(Dof*(i-1)+1:Dof*i))
          do j=1,DoF
             x=(i-1)*dx+(j-1)*ddx
             write(2,*),x,vector_work(j)
          end do
       end do
    else
       do i=1,nb_elem
          do j=1,DoF
             x=(i-1)*dx+(j-1)*ddx
             write(2,*),x,vector((i-1)*DoF+j)
          end do
       end do
    end if
  end subroutine print_vect


  subroutine print_sol(problem,N)
    type(acoustic_problem),intent(in) :: problem
    integer               ,intent(in) :: N

    call print_vect(problem%U,problem%nb_elem,problem%DoF,problem%dx,           &
                    problem%bernstein,N,'U')

    call print_vect(problem%P,problem%nb_elem,problem%DoF,problem%dx,           &
                    problem%bernstein,N,'P')

  end subroutine print_sol

  
end module m_acoustic
