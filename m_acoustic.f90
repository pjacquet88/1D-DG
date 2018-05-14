module m_acoustic
  use m_polynom
  use m_matrix
  use m_powermethod
  implicit none

  type acoustic_problem
     integer                         :: DoF
     integer                         :: nb_elem
     logical                         :: bernstein
     logical                         :: F_forte
     real                            :: dx
     real                            :: dt
     real                            :: total_length
     real                            :: final_time
     real                            :: alpha
     real,dimension(:),allocatable   :: U
     real,dimension(:),allocatable   :: P
     character(len=20)               :: boundaries
     character(len=20)               :: signal
     type(sparse_matrix)             :: Ap
     type(sparse_matrix)             :: Av
     real,dimension(:,:),allocatable :: Minv
     integer                         :: k_max
     real                            :: epsilon
  end type acoustic_problem

  real,dimension(15)              :: xx,weight

  real,parameter :: PI=acos(-1.0)

  public  :: init_problem,one_time_step,init_ApAv,                              &
             solution,free_problem,                                             &
             print_vect,print_sol,error_periodique

  private :: xx,weight,                                                         &
             init_quadrature,Int_bxb,Int_lxl,init_m_loc,init_s_loc,             &
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
    real,dimension(:),intent(in) :: bpol1
    real,dimension(:),intent(in) :: bpol2
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
    real,dimension(:),intent(in) :: lpol1
    real,dimension(:),intent(in) :: lpol2
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
  
  function solution(x,dt,a)
    real,intent(in) :: x
    real,intent(in) :: dt
    character(len=*),intent(in) :: a
    real                        :: solution

    real                        :: xt

    xt=x+dt
    
    if (a.eq.'sinus') then
       solution=sin(4*PI*xt)
       ! solution=sin(xt)
    else if (a.eq.'creneau') then
       if ((x.lt.0.55).and.(x.gt.0.45)) then
          solution=0.5
       else
          solution=0.0
       end if
    else if (a.eq.'x') then
       solution=xt**5
    else
       solution=(xt-0.5)*exp(-(2.0*PI*(xt-0.5)*2.0)**2.0)*20
    end if    
  end function solution


  subroutine init_problem(problem,nb_elem,DoF,total_length,final_time,alpha,    &
                          bernstein,F_forte,signal,boundaries,k_max,epsilon)
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
    integer               ,intent(in)  :: k_max
    real                  ,intent(in)  :: epsilon
    
    integer :: i,j
    real    :: gd,dd

    problem%DoF=DoF
    problem%nb_elem=nb_elem
    problem%total_length=total_length
    problem%dx=total_length/(nb_elem)
    problem%alpha=alpha
    problem%bernstein=bernstein
    problem%F_forte=F_forte
    problem%boundaries=boundaries
    problem%k_max=k_max
    problem%epsilon=epsilon
    problem%signal=signal

    allocate(problem%U(DoF*nb_elem),problem%P(DoF*nb_elem))
    problem%U=0.0
    problem%P=0.0

    print*,'nb_elem              ::',nb_elem
    print*,'DoF                  ::',DoF
    print*,'total_length         ::',total_length
    print*,'final_time           ::',final_time
    print*,'dx                   ::',problem%dx
    print*,'dt                   ::',problem%dt
    print*,'Boundary Condition   ::','   ',boundaries
    print*,'Penalisation alpha   ::',alpha
    print*,'Max Iter Power-Method::',k_max
    print*,'Epsilon Power-Method ::',epsilon
    print*,' '
    if (bernstein) then
       print*,'The problem is solved with Bernstein elements'
    else
       print*,'The problem is solved with Lagrange elements'
    end if
    if (F_forte) then
       print*,'The problem is solved with a strong formulation'
    else
       print*,'The problem is solved with a weak formulation'
    end if
    print*,'------------------------------------------------------------'    
  end subroutine init_problem

  
  subroutine free_problem(problem)
    type(acoustic_problem),intent(inout) :: problem
    deallocate(problem%U)
    deallocate(problem%P)
    deallocate(problem%Minv)
  end subroutine free_problem

  
  subroutine init_UP(problem)
    type(acoustic_problem),intent(inout) :: problem
    integer                              :: i,j,DoF
    real    :: x,ddx
    
    DoF=problem%DoF
    ddx=problem%dx/(DoF-1)
    x=0.0
    do i=1,problem%nb_elem
       do j=1,DoF
          x=(i-1)*problem%dx+(j-1)*ddx
          problem%U((i-1)*problem%DoF+j)=solution(x,0.0,problem%signal)
          problem%P((i-1)*problem%DoF+j)=solution(x,-problem%dt/2,problem%signal)
       end do
    end do

    if (problem%bernstein) then
       do i=1,problem%nb_elem
          problem%U(DoF*(i-1)+1:Dof*i)=matmul(L2B,problem%U(DoF*(i-1)+1:DoF*i))
          problem%P(DoF*(i-1)+1:Dof*i)=matmul(L2B,problem%P(DoF*(i-1)+1:DoF*i))
       end do
    end if
  end subroutine init_UP
  

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
          m_loc(i,i)=Int_bxb(base_b(i,:),base_b(i,:))
          do j=i+1,DoF
             m_loc(i,j)=Int_bxb(base_b(i,:),base_b(j,:))
             m_loc(j,i)=m_loc(i,j)
          end do
       end do
    else
       do i=1,DoF
          m_loc(i,i)=Int_lxl(base_l(i,:),base_l(i,:))
          do j=i+1,DoF
             m_loc(i,j)=Int_lxl(base_l(i,:),base_l(j,:))
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
    real,dimension(:),allocatable :: bpol,dbpol
    real,dimension(:),allocatable :: dlpol

    DoF=size(s_loc,1)
    
    allocate(bpol(DoF))       
    allocate(dbpol(DoF))       
    allocate(dlpol(DoF))          
    
    if (bernstein) then
       s_loc=D1-D0
    else
       do i=1,DoF
          do j=1,DoF
             call Lagrange2Bernstein(base_l(j,:),bpol)
             call deriv_pol_b(bpol,dbpol)
             call Bernstein2Lagrange(dbpol,dlpol,1,DoF)
             s_loc(i,j)=Int_lxl(base_l(i,:),dlpol)
          end do
       end do
    end if
    deallocate(bpol)      
    deallocate(dbpol)
    deallocate(dlpol) 
  end subroutine init_s_loc


  subroutine init_dt(Ap,Av,dt,k_max,epsilon,nb_elem)
    real,dimension(:,:),intent(in)  :: Ap
    real,dimension(:,:),intent(in)  :: Av
    real               ,intent(out) :: dt
    integer            ,intent(in)  :: k_max
    real               ,intent(in)  :: epsilon
    integer            ,intent(in)  :: nb_elem
    type(sparse_matrix)             :: sparse_A
    real,dimension(size(Ap,1),size(Ap,2)) :: A
    real                            :: max_value
    integer                         :: i
    real,parameter                  :: alpha=1.90
    real                            :: cfl

    ! A=matmul(Ap,Av)
    ! ! do i=1,size(A,1)
    ! !    write(789,*) A(i,:)
    ! ! end do
    
    ! call Full2Sparse(A,sparse_A)
    ! call power_method_sparse(sparse_A,max_value,k_max,epsilon)
    ! dt=alpha/(sqrt(abs(max_value)))
    ! call free_sparse_matrix(sparse_A)
    ! ! A=matmul(Av,Ap)
    ! ! call Full2Sparse(matmul(Av,Ap),sparse_A)
    ! ! call power_method_sparse(sparse_A,max_value,k_max,epsilon)
    ! ! dt=min(dt,alpha/(sqrt(abs(max_value))))
    ! ! call free_sparse_matrix(sparse_A)



    !****** dt constant *************************
    !****** 100 ***************************
    ! dt=4.16666589E-04*alpha  !ordre  1 
    ! dt= 2.08152022E-04*alpha  !ordre  2
    ! dt= 1.24648141E-04*alpha  !ordre  3
     dt= 8.29414494E-05*alpha  !ordre  4
    ! dt= 5.91463395E-05*alpha  !ordre  5
    ! dt= 4.42988749E-05*alpha  !ordre  6
    ! dt= 3.44150467E-05*alpha  !ordre  7
    ! dt= 2.75054135E-05*alpha  !ordre  8
    ! dt= 2.24858413E-05*alpha  !ordre  9
    ! dt= 1.87254609E-05*alpha  !ordre  10

    !****** 500 ****************************
    ! dt= 1.65833015E-04        !ordre  1 
    ! dt= 8.28442862E-05        !ordre  2
    ! dt= 4.96104076E-05        !ordre  3
    ! dt= 8.29414494E-05*alpha  !ordre  4
    ! dt= 5.91463395E-05*alpha  !ordre  5
    ! dt= 4.42988749E-05*alpha  !ordre  6
    ! dt= 3.44150467E-05*alpha  !ordre  7
    ! dt= 2.75054135E-05*alpha  !ordre  8
    ! dt= 2.24858413E-05*alpha  !ordre  9
    ! dt= 1.87254609E-05*alpha  !ordre  10
     !*********************************************

     !************** cfl constante ****************
     cfl=0.95/(dt*100)
     dt=1.0/(cfl*nb_elem)

    
  end subroutine init_dt
  

  subroutine init_ApAv(problem)
    type(acoustic_problem),intent(inout)        :: problem

    real,dimension(problem%nb_elem*problem%DoF,                                 &
                   problem%nb_elem*problem%DoF) :: s_glob,MINV,F1,F2,           &
                                                   Av_full,Ap_full
    real,dimension(problem%DoF,problem%DoF)     :: m_loc
    real,dimension(problem%DoF,problem%DoF)     :: s_loc
    real                                        :: alpha
    
    integer                 :: i,j
    integer                 :: DoF,nb_elem

    DoF=problem%DoF           ! For sake of lisibility
    nb_elem=problem%nb_elem

    Ap_full=0.0
    Av_full=0.0
    F1=0.0
    F2=0.0

    alpha=problem%alpha
    
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

       alpha=-problem%alpha
       
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
           
    if (problem%boundaries.eq.'periodique') then
       alpha=problem%alpha

       F1(1,1)=0.5+alpha
       F1(1,DoF*nb_elem)=-0.5-alpha

       F1(Dof*nb_elem,Dof*nb_elem)=-0.5+alpha
       F1(Dof*nb_elem,1)=0.5-alpha

       alpha=-problem%alpha

       F2(1,1)=0.5+alpha
       F2(1,DoF*nb_elem)=-0.5-alpha

       F2(Dof*nb_elem,Dof*nb_elem)=-0.5+alpha
       F2(Dof*nb_elem,1)=0.5-alpha

    elseif (problem%boundaries.eq.'dirichlet') then
       alpha=problem%alpha

       F1(1,1)=1.0+2.0*alpha
       F1(1,DoF*nb_elem)=0.0

       F1(Dof*nb_elem,Dof*nb_elem)=-1.0-2.0*alpha
       F1(Dof*nb_elem,1)=0.0

       alpha=-problem%alpha

       F2(1,1)=0.0
       F2(1,DoF*nb_elem)=0.0

       F2(Dof*nb_elem,Dof*nb_elem)=0.0
       F2(Dof*nb_elem,1)=0.0
    elseif (problem%boundaries.eq.'neumann') then
     alpha=problem%alpha

       F1(1,1)=0.0
       F1(1,DoF*nb_elem)=0.0

       F1(Dof*nb_elem,Dof*nb_elem)=0.0
       F1(Dof*nb_elem,1)=0.0

       alpha=-problem%alpha

       F2(1,1)=1.0-2.0*alpha
       F2(1,DoF*nb_elem)=0.0

       F2(Dof*nb_elem,Dof*nb_elem)=-1.0+2.0*alpha
       F2(Dof*nb_elem,1)=0.0
    end if

    call init_quadrature
    allocate(problem%Minv(DoF,DoF))
    call init_m_loc(m_loc,problem%Minv,problem%bernstein)
    call init_s_loc(s_loc,problem%bernstein)

    MINV=0.0
    s_glob=0.0
    do i=1,nb_elem
       MINV(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=problem%Minv
    end do
    do i=1,nb_elem
       s_glob(DoF*(i-1)+1:Dof*i,DoF*(i-1)+1:Dof*i)=s_loc
    end do

    if (problem%bernstein) then
       Ap_full=(s_glob+matmul(MINV,F1))
       Av_full=(s_glob+matmul(MINV,F2))

       if (.not.problem%F_forte) then
          Ap_full=-transpose(Ap_full)
          Av_full=-transpose(Av_full)
       end if

       call Full2Sparse(Ap_full,problem%Ap)
       call Full2Sparse(Av_full,problem%Av)

       print*,'---------------------------'
       print*,'epsilon',problem%epsilon
       print*,'k_max',problem%k_max
       print*,'---------------------------'

       call init_dt((1/problem%dx)*Ap_full,(1/problem%dx)*Av_full,problem%dt,   &
                     problem%k_max,problem%epsilon,problem%nb_elem)
       
       problem%Ap%Values=(problem%dt/problem%dx)*problem%Ap%Values
       problem%Av%Values=(problem%dt/problem%dx)*problem%Av%Values
    else
          Ap_full=s_glob+F1
          Av_full=s_glob+F2   

       if (.not.problem%F_forte) then
          Ap_full=-transpose(Ap_full)
          Av_full=-transpose(Av_full)
       end if

       Ap_full=matmul(MINV,Ap_full)
       Av_full=matmul(MINV,Av_full)

       call Full2Sparse(Ap_full,problem%Ap)
       call Full2Sparse(Av_full,problem%Av)
       
       call init_dt((1/problem%dx)*Ap_full,(1/problem%dx)*Av_full,problem%dt,   &
                     problem%k_max,problem%epsilon,problem%nb_elem)
       
       problem%Ap%Values=(problem%dt/problem%dx)*problem%Ap%Values
       problem%Av%Values=(problem%dt/problem%dx)*problem%Av%Values
    end if
    
    print*,'Non-Zero Values of Ap,Av  ::',problem%Ap%NNN
    print*,'Size of Ap,AV matrices    ::',problem%Ap%nb_ligne,'x', &
         problem%Ap%nb_ligne,'=',problem%Ap%nb_ligne**2
    print*,'Ratio                     ::',real(problem%Ap%NNN)/problem%Ap%nb_ligne**2
  end subroutine init_ApAv
  

  subroutine one_time_step(problem,t)
    type(acoustic_problem),intent(inout) :: problem
    real,                  intent(in)    :: t

    real,dimension(problem%DoF*problem%nb_elem) :: RHSv,RHSp
    real                                        :: gd,gg

    RHSp=0.0
    RHSv=0.0

    if (problem%boundaries.eq.'dirichlet') then
       gg=min(0.5*t,0.5)
       gd=-gg
       RHSp(1)=gg*(-1.0-2.0*problem%alpha)*(problem%dt/problem%dx)
       RHSp(size(RHSp))=gd*(1.0+2.0*problem%alpha)*(problem%dt/problem%dx)
       
       RHSv(1)=0.0
       RHSv(size(RHSv))=0.0

       RHSp(1:problem%DoF)=matmul(problem%Minv,RHSp(1:problem%DoF))
       RHSp(size(RHSp)-problem%DoF+1:size(RHSp))=matmul(problem%Minv,           &
            RHSp(size(RHSp)-problem%DoF+1:size(RHSp)))
    end if

       problem%U=problem%U-sparse_matmul(problem%Av,problem%P)-RHSv
       problem%P=problem%P-sparse_matmul(problem%Ap,problem%U)-RHSp  
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
             write(2,*)x,vector_work(j)
          end do
       end do
    else
       do i=1,nb_elem
          do j=1,DoF
             x=(i-1)*dx+(j-1)*ddx
             write(2,*)x,vector((i-1)*DoF+j)
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

  subroutine error_periodique(problem,t,errorU,errorP,iter)
    type(acoustic_problem),intent(in)  :: problem
    real                  ,intent(in)  :: t
    real                  ,intent(out) :: errorU
    real                  ,intent(out) :: errorP
    integer,optional      ,intent(in)  :: iter
    
    real,dimension(problem%nb_elem*problem%DoF) :: U_ex
    real,dimension(problem%nb_elem*problem%DoF) :: U_work
    real,dimension(problem%nb_elem*problem%DoF) :: P_ex
    real,dimension(problem%nb_elem*problem%DoF) :: P_work
    real                                        :: dx,ddx,x
    integer                                     :: i,j,k
    integer                                     :: DoF,nb_elem

    DoF=problem%DoF
    nb_elem=problem%nb_elem

    dx=problem%total_length/nb_elem
    ddx=dx/(DoF-1)
    
    do i=1,nb_elem
       do j=1,DoF
          x=(i-1)*problem%dx+(j-1)*ddx-t

          do while (x.le.0.0)
             x=x+problem%total_length
          end do
          
          U_ex((i-1)*DoF+j)=solution(x,0.0,problem%signal)
          P_ex((i-1)*DoF+j)=solution(x,-problem%dt/2,problem%signal)
       end do
    end do

    errorU=0
    errorP=0

    if(problem%bernstein) then
       do k=1,nb_elem
          !U_work(DoF*(k-1)+1:DoF*k)=matmul(B2L,problem%U(DoF*(k-1)+1:DoF*k))
          !P_work(DoF*(k-1)+1:DoF*k)=matmul(B2L,problem%P(DoF*(k-1)+1:DoF*k))
          call Bernstein2Lagrange(problem%U,U_work,nb_elem,DoF)
          call Bernstein2Lagrange(problem%P,P_work,nb_elem,DoF)
       end do
    else
       U_work=problem%U
       P_work=problem%P
    end if
    do i=1,nb_elem*DoF
       errorU=errorU+(U_ex(i)-U_work(i))**2
       errorP=errorP+(P_ex(i)-P_work(i))**2
    end do
    if (present(iter)) then
       call print_vect(U_ex-U_work,nb_elem,DoF,problem%dx,.false.,iter,'Error')
    else
       call print_vect(0.0*U_work,0,DoF,problem%dx,.false.,0,'Error')
       call print_vect(0.0*U_work,1,DoF,problem%dx,.false.,1,'Error')
       call print_vect(U_ex-U_work,nb_elem,DoF,problem%dx,.false.,nb_elem,'Error')
    end if
    errorU=sqrt(errorU)/sqrt(dot_product(U_ex,U_ex))
    errorP=sqrt(errorP)/sqrt(dot_product(P_ex,P_ex))

    print*,'-----------------------------'
    print*,'-----------------------------'
    print*,'t=',t
    print*,'Lerreur en U est de : ', errorU
    print*,'Lerreur en P est de : ', errorP
    print*,'-----------------------------'    
    print*,'-----------------------------'
    
    if (present(iter)) then
       call print_vect(U_ex,nb_elem,DoF,problem%dx,           &
            .false.,iter,'Uex')
       ! call print_vect(P_ex,nb_elem,DoF,problem%dx,           &
       !      .false.,0,'Pex')
    end if
  end subroutine error_periodique
  
end module m_acoustic
