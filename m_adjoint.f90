module m_adjoint
  use m_polynom
  use m_matrix
  use m_powermethod
  use m_file_function
  use m_time_scheme
  implicit none

  ! Acoustinc Adjoint problem type
  type adjoint_problem
     integer                         :: DoF          ! nb of degree of freedom 
     integer                         :: nb_elem      ! nb of elements
     character(len=20)               :: time_scheme  ! the used time scheme
     logical                         :: bernstein    ! bool for bernstein elmnts
     real,dimension(:) ,allocatable  :: density      ! density model
     real,dimension(:) ,allocatable  :: velocity     ! velocity model
     real                            :: dx           ! element length
     real                            :: dt           ! step time length
     real                            :: total_length ! domain length
     real                            :: final_time   ! final time
     real                            :: alpha        ! penalisation value
     real,dimension(:),allocatable   :: U            ! velocity unknowns vector
     real,dimension(:),allocatable   :: P            ! pressure unknowns vector
     real,dimension(:),allocatable   :: Pk1,Pk2      ! AB3 vectors
     real,dimension(:),allocatable   :: Uk1,Uk2      ! AB3 vectors
     character(len=20)               :: boundaries   ! string for BC
     character(len=20)               :: signal       ! string to initialize U,P
     !----------- MATRICES ---------------------------
     type(sparse_matrix)             :: Ap,App       
     type(sparse_matrix)             :: Av
     type(sparse_matrix)             :: Minv_p,Minv_v
     type(sparse_matrix)             :: M,Minv
     !----------- SOURCE & RECEIVER ------------------
     integer                         :: source_loc    ! beginning of an element
     integer                         :: receiver_loc  ! beginning of an element
     integer                         :: n_time_step   ! nb of time step

     real,dimension(:,:),allocatable :: data_P        ! pressure data received through time
     real,dimension(:,:),allocatable :: data_U        ! velocity data received through time
     real,dimension(:,:),allocatable :: P_received    ! simulated pressure received
     real,dimension(:,:),allocatable :: U_received    ! simulated data received

     character(len=20)               :: strategy      ! Strategy employed : ATD or DTA 
     character(len=20)               :: scalar_product! canonical or M inner product

     real,dimension(:)  ,allocatable :: RHSu
     real,dimension(:)  ,allocatable :: RHSu_half
     real,dimension(:)  ,allocatable :: RHSp
     real,dimension(:)  ,allocatable :: RHSp_half
  end type adjoint_problem
  
  real,dimension(15)                 :: xx,weight
  real,parameter :: PI=acos(-1.0)

  public  :: init_adjoint_problem,init_adjoint_operator,all_backward_step,      &
             one_backward_step,free_adjoint_problem,Int_lxl,Int_bxb

  private :: xx,weight,initial_perturbation,PI,                                             &
             init_quadrature,init_m_loc,init_s_loc,                             &
             init_UP,init_m_glob,init_minv_glob,init_minv_glob_abc,init_mabs,   &
             init_s_glob,init_FpFv,init_DpDv,init_dt,eval_backward_source
  
contains


  !****************** QUADRATURE ***********************************************
  ! Initializes the quadrature to evaluate polynomial integrals
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

  ! Evaluates the integral of the product of 2 bernstein polyniomials
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


  ! Evaluates the integral of the product of 2 lagrange polynomials
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
  ! Contains the initial perturbation
  function initial_perturbation(x,dt,signal)
    real,intent(in) :: x
    real,intent(in) :: dt
    character(len=*),intent(in) :: signal
    real                        :: initial_perturbation

    real                        :: xt

    xt=x+dt
    
    if (signal.eq.'sinus') then
       initial_perturbation=sin(4*PI*xt)
       ! initial_perturbation=sin(xt)
    else if (signal.eq.'boxcar') then
       if ((x.lt.0.55).and.(x.gt.0.45)) then
          initial_perturbation=0.5
       else
          initial_perturbation=0.0
       end if
    else if (signal.eq.'x') then
       initial_perturbation=xt**5
        else if (signal.eq.'flat') then
       initial_perturbation=0.0
    else
       initial_perturbation=(xt-0.5)*exp(-(2.0*PI*(xt-0.5)*2.0)**2.0)*20
    end if    
  end function initial_perturbation


  ! Initializes the adjoint acoustic problem
  subroutine init_adjoint_problem(problem,nb_elem,DoF,time_scheme,velocity,     &
                                  density,total_length,final_time,              &
                                  alpha,bernstein,signal,boundaries,source_loc, &
                                  receiver_loc, data_P,data_U,P_received,       &
                                  U_received,strategy,scalar_product)
    
    type(adjoint_problem) ,intent(out) :: problem
    integer               ,intent(in)  :: nb_elem
    integer               ,intent(in)  :: DoF
    real,dimension(:)     ,intent(in)  :: velocity
    real,dimension(:)     ,intent(in)  :: density
    character(len=*)      ,intent(in)  :: time_scheme
    real                  ,intent(in)  :: total_length
    real                  ,intent(in)  :: final_time
    real                  ,intent(in)  :: alpha
    logical               ,intent(in)  :: bernstein
    character(len=*)      ,intent(in)  :: signal
    character(len=*)      ,intent(in)  :: boundaries
    integer               ,intent(in)  :: source_loc
    integer               ,intent(in)  :: receiver_loc
    real,dimension(:,:)   ,intent(in)  :: data_P
    real,dimension(:,:)   ,intent(in)  :: data_U
    real,dimension(:,:)   ,intent(in)  :: P_received
    real,dimension(:,:)   ,intent(in)  :: U_received
    character(len=*)      ,intent(in)  :: strategy
    character(len=*)      ,intent(in)  :: scalar_product
    
    integer :: i,j
    integer :: dv,dp,size_velocity,size_density
    real    :: gd,dd

    problem%DoF=DoF
    problem%time_scheme=time_scheme
    problem%nb_elem=nb_elem
    problem%total_length=total_length
    problem%dx=total_length/(nb_elem)
    problem%final_time=final_time
    problem%alpha=alpha
    problem%bernstein=bernstein
    problem%boundaries=boundaries
    problem%signal=signal
    problem%receiver_loc=receiver_loc
    problem%source_loc=source_loc
    
    size_velocity=size(velocity)
    size_density=size(density)
    
    allocate(problem%U(DoF*nb_elem),problem%P(DoF*nb_elem))
    problem%U=0.0
    problem%P=0.0
    
    if (problem%time_scheme.eq.'AB3') then
       allocate(problem%Pk1(DoF*nb_elem),problem%Uk1(DoF*nb_elem))
       allocate(problem%Pk2(DoF*nb_elem),problem%Uk2(DoF*nb_elem))
       problem%Pk1=0.0
       problem%Pk2=0.0
       problem%Uk1=0.0
       problem%Uk2=0.0
    end if
    
    if (size(velocity).gt.nb_elem) then
       print*,'Velocity input size does not correspond'
       STOP
    end if
    if (size(density).gt.nb_elem) then
       print*,'Density input size does not correspond'
       STOP
    end if
    allocate(problem%velocity(nb_elem),problem%density(nb_elem))

    dv=max(nb_elem/size(velocity),1)
    do i=1,size(velocity)
       do j=(i-1)*dv+1,i*dv
          problem%velocity(j)=velocity(i)
       end do
    end do
    do j=size(velocity)*dv+1,nb_elem
       problem%velocity(j)=velocity(size(velocity))
    end do
    
    dp=max(nb_elem/size(density),1)
    do i=1,size(density)
       do j=(i-1)*dp+1,i*dp
          problem%density(j)=density(i)
       end do
    end do
    do j=size(density)*dp+1,nb_elem
       problem%density(j)=density(size(density))
    end do

    problem%strategy=strategy
    problem%scalar_product=scalar_product
    
    
    !---------------init matrices------------------
    if (problem%boundaries.eq.'ABC') then
       call init_adjoint_operator_abc(problem)
    else
       call init_adjoint_operator(problem)
    end if
    
    problem%n_time_step=int(problem%final_time/problem%dt)+1
    problem%final_time=problem%dt*problem%n_time_step
    ! here I change the final time to fit frames of the forward and the backward
    
    allocate(problem%data_P(size(data_P,1),size(data_P,2)))
    problem%data_P=data_P
    allocate(problem%data_U(size(data_U,1),size(data_U,2)))
    problem%data_U=data_U
    allocate(problem%P_received(size(P_received,1),size(P_received,2)))
    problem%P_received=P_received
    allocate(problem%U_received(size(U_received,1),size(U_received,2)))   
    problem%U_received=U_received

    allocate(problem%RHSu(nb_elem*DoF))
    allocate(problem%RHSp(nb_elem*DoF))
    if (time_scheme.eq.'RK4') then
       allocate(problem%RHSu_half(nb_elem*DoF))
       allocate(problem%RHSp_half(nb_elem*DoF))
    end if
  end subroutine init_adjoint_problem


  ! Initialize P and U at time=0
  subroutine init_UP(problem)
    type(adjoint_problem),intent(inout) :: problem
    integer                              :: i,j,DoF
    real    :: x,ddx,dephasing
    
    DoF=problem%DoF
    ddx=problem%dx/(DoF-1)
    x=0.0

    if ((problem%time_scheme.eq.'LF').or.(problem%time_scheme.eq.'LF4')) then
       dephasing=-problem%dt/2.0
    else
       dephasing=0.0
    end if

    do i=1,problem%nb_elem
       do j=1,DoF
          x=(i-1)*problem%dx+(j-1)*ddx
          problem%U((i-1)*problem%DoF+j)=initial_perturbation(x,0.0,problem%signal)
          problem%P((i-1)*problem%DoF+j)=initial_perturbation(x,dephasing,problem%signal)
       end do
    end do

    if (problem%bernstein) then
       do i=1,problem%nb_elem
          problem%U(DoF*(i-1)+1:Dof*i)=matmul(L2B,problem%U(DoF*(i-1)+1:DoF*i))
          problem%P(DoF*(i-1)+1:Dof*i)=matmul(L2B,problem%P(DoF*(i-1)+1:DoF*i))
       end do
    end if
  end subroutine init_UP
  

 !--------------------- Init Matrices ------------------------------------
  ! Initializes the local masse matrix (on one element)
  subroutine init_m_loc(m_loc,DoF,bernstein)
    real,dimension(:,:),intent(out) :: m_loc
    integer            ,intent(in)  :: DoF
    logical            ,intent(in)  :: bernstein
    integer                         :: i,j

    m_loc=0.0
    
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
    
  end subroutine init_m_loc

  ! Initializes the local inverse mass matrix (on one element)
  subroutine init_minv_loc(minv_loc,m_loc)
    real,dimension(:,:),intent(out) :: minv_loc
    real,dimension(:,:),intent(in)  :: m_loc
    minv_loc=0.0
    minv_loc=LU_inv(m_loc)
  end subroutine init_minv_loc

  ! Initialize the global mass matrix (on all the elements)  
  subroutine init_m_glob(m_glob,m_loc,nb_elem)
    real,dimension(:,:),intent(out) :: m_glob
    real,dimension(:,:),intent(in)  :: m_loc
    integer            ,intent(in)  :: nb_elem

    integer :: i,DoF

    m_glob=0.0
    
    DoF=size(m_loc,1)
        
    do i=1,nb_elem
       m_glob(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=m_loc
    end do
  end subroutine init_m_glob
  

  ! Initialize the global inverse mass matrix (on all the element)
  subroutine init_minv_glob(minv_glob,minv_loc,nb_elem)
    real,dimension(:,:),intent(out) :: minv_glob       
    real,dimension(:,:),intent(in)  :: minv_loc
    integer            ,intent(in)  :: nb_elem
    integer                         :: i,DoF

    minv_glob=0.0

    DoF=size(minv_loc,1)
    
    do i=1,nb_elem
       minv_glob(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=minv_loc
    end do

  end subroutine init_minv_glob


  ! initialize the global inverse "abc mass matrix" 
  subroutine init_minv_glob_abc(minv_glob_abc,m_loc,nb_elem,DoF,velocity_beg,   &
                                velocity_end,dx,dt,time_scheme)
    real,dimension(:,:),intent(out) :: minv_glob_abc
    real,dimension(:,:),intent(in)  :: m_loc
    integer            ,intent(in)  :: nb_elem
    integer            ,intent(in)  :: DoF
    real               ,intent(in)  :: velocity_beg
    real               ,intent(in)  :: velocity_end
    real               ,intent(in)  :: dx
    real               ,intent(in)  :: dt
    character(len=*)   ,intent(in)  :: time_scheme

    integer                         :: i
    real,dimension(DoF,DoF)         :: m_loc_abc1,m_loc_abc2
    real,dimension(DoF,DoF)         :: minv_loc,minv_loc_abc1,minv_loc_abc2

    minv_glob_abc=0.0

    m_loc_abc1=m_loc
    m_loc_abc2=m_loc

    if (time_scheme.eq.'LF') then
       m_loc_abc2(DoF,DoF)=m_loc_abc2(DoF,DoF)+0.5*velocity_end*dt/dx
       m_loc_abc1(1,1)=m_loc_abc1(1,1)+0.5*velocity_beg*dt/dx
    end if

    minv_loc=LU_inv(m_loc)
    minv_loc_abc1=LU_inv(m_loc_abc1)
    minv_loc_abc2=LU_inv(m_loc_abc2)
    
    minv_glob_abc=0.0
    minv_glob_abc(1:DoF,1:DoF)=minv_loc_abc1
    do i=2,nb_elem-1
       minv_glob_abc(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=minv_loc
    end do
    minv_glob_abc(DoF*(nb_elem-1)+1:DoF*nb_elem,DoF*(nb_elem-1)+1:DoF*nb_elem)  &
         =minv_loc_abc2
  end subroutine init_minv_glob_abc


  ! Initializes the global Abs matrix
  subroutine init_mabs(mabs,nb_elem,DoF,velocity_beg,velocity_end,dx,dt,time_scheme)
    real,dimension(:,:),intent(out) :: mabs
    integer            ,intent(in)  :: nb_elem
    integer            ,intent(in)  :: DoF
    real               ,intent(in)  :: velocity_beg
    real               ,intent(in)  :: velocity_end
    real               ,intent(in)  :: dx
    real               ,intent(in)  :: dt
    character(len=*)   ,intent(in)  :: time_scheme

    integer :: i

    mabs=0.0

    if (time_scheme.eq.'LF') then
       mabs(DoF*nb_elem,DoF*nb_elem)=0.5*velocity_end*dt/dx
       mabs(1,1)=0.5*velocity_beg*dt/dx
    else
       mabs(DoF*nb_elem,DoF*nb_elem)=velocity_end*dt/dx
       mabs(1,1)=velocity_beg*dt/dx
    end if
  end subroutine init_mabs


  ! Initailizes the local stiffness matrix (on one element)
  subroutine init_s_loc(s_loc,DoF,bernstein)
    real,dimension(:,:),intent(out) :: s_loc
    integer            ,intent(in)  :: DoF
    logical            ,intent(in)  :: bernstein
    integer                         :: i,j
    real,dimension(DoF)             :: bpol,dbpol
    real,dimension(DoF)             :: dlpol

    s_loc=0.0
    
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
  end subroutine init_s_loc


  ! Initialize the global stiffness matrix (on all the elements)
  subroutine init_s_glob(s_glob,s_loc,nb_elem)
    real,dimension(:,:),intent(out) :: s_glob
    real,dimension(:,:),intent(in)  :: s_loc
    integer            ,intent(in)  :: nb_elem
    integer :: i,DoF

    DoF=size(s_loc,1)
    s_glob=0.0
    
    do i=1,nb_elem
       s_glob(DoF*(i-1)+1:Dof*i,DoF*(i-1)+1:Dof*i)=s_loc
    end do
  end subroutine init_s_glob

  ! Initalize the Fp and Fv flux matrices (with penalization)
  subroutine init_FpFv(Fp,Fv,DoF,nb_elem,boundaries,alpha,receiver_loc)
    real,dimension(:,:),intent(out) :: Fp
    real,dimension(:,:),intent(out) :: Fv
    integer            ,intent(in)  :: DoF
    integer            ,intent(in)  :: nb_elem
    character(len=*)   ,intent(in)  :: boundaries
    real               ,intent(in)  :: alpha
    integer            ,intent(in)  :: receiver_loc

    integer :: i,j
    real    :: alpha_dummy

    Fp=0.0
    Fv=0.0

    alpha_dummy=alpha

    Fp(DoF,DoF)=-0.5+alpha_dummy
    Fp(DoF,DoF+1)=0.5-alpha_dummy

    do i=2,nb_elem-1
       Fp(Dof*(i-1)+1,Dof*(i-1)+1)=0.5+alpha_dummy
       Fp(Dof*(i-1)+1,Dof*(i-1))=-0.5-alpha_dummy

       Fp(Dof*i,Dof*i)=-0.5+alpha_dummy
       Fp(Dof*i,Dof*i+1)=0.5-alpha_dummy
    end do
    Fp(Dof*(nb_elem-1)+1,Dof*(nb_elem-1)+1)=0.5+alpha_dummy
    Fp(Dof*(nb_elem-1)+1,Dof*(nb_elem-1))=-0.5-alpha_dummy

    alpha_dummy=-alpha

    Fv(DoF,DoF)=-0.5+alpha_dummy
    Fv(DoF,DoF+1)=0.5-alpha_dummy

    do i=2,nb_elem-1
       Fv(Dof*(i-1)+1,Dof*(i-1)+1)=0.5+alpha_dummy
       Fv(Dof*(i-1)+1,Dof*(i-1))=-0.5-alpha_dummy

       Fv(Dof*i,Dof*i)=-0.5+alpha_dummy
       Fv(Dof*i,Dof*i+1)=0.5-alpha_dummy
    end do
    Fv(Dof*(nb_elem-1)+1,Dof*(nb_elem-1)+1)=0.5+alpha_dummy
    Fv(Dof*(nb_elem-1)+1,Dof*(nb_elem-1))=-0.5-alpha_dummy

    if (boundaries.eq.'periodic') then
       alpha_dummy=alpha

       Fp(1,1)=0.5+alpha_dummy
       Fp(1,DoF*nb_elem)=-0.5-alpha_dummy

       Fp(Dof*nb_elem,Dof*nb_elem)=-0.5+alpha_dummy
       Fp(Dof*nb_elem,1)=0.5-alpha_dummy
       
       alpha_dummy=-alpha

       Fv(1,1)=0.5+alpha_dummy
       Fv(1,DoF*nb_elem)=-0.5-alpha_dummy

       Fv(Dof*nb_elem,Dof*nb_elem)=-0.5+alpha_dummy
       Fv(Dof*nb_elem,1)=0.5-alpha_dummy

    elseif (boundaries.eq.'dirichlet') then

       alpha_dummy=alpha

       Fp(1,1)=0.0
       Fp(1,DoF*nb_elem)=0.0

       Fp(Dof*nb_elem,Dof*nb_elem)=0.0
       Fp(Dof*nb_elem,1)=0.0


       alpha_dummy=-alpha

       Fv(1,1)=1.0!+2.0*alpha_dummy
       Fv(1,DoF*nb_elem)=0.0

       Fv(Dof*nb_elem,Dof*nb_elem)=-1.0!-2.0*alpha_dummy
       Fv(Dof*nb_elem,1)=0.0


    elseif (boundaries.eq.'ABC') then

       alpha_dummy=alpha

       Fv(1,1)=0.0!(1.0+2.0*alpha_dummy)
       Fv(1,DoF*nb_elem)=0.0

       Fv(Dof*nb_elem,Dof*nb_elem)=0.0!-1.0-2.0*alpha_dummy
       Fv(Dof*nb_elem,1)=0.0
       
       alpha_dummy=-alpha

!       Fp(1,1)=1.0
       Fp(1,1)=1.0
       Fp(1,DoF*nb_elem)=0.0

       Fp(Dof*nb_elem,Dof*nb_elem)=-1.0!*2
       Fp(Dof*nb_elem,1)=0.0

    elseif (boundaries.eq.'neumann') then
       alpha_dummy=alpha
       
       Fp(1,1)=1.0!-2.0*alpha_dummy
       Fp(1,DoF*nb_elem)=0.0
       
       Fp(Dof*nb_elem,Dof*nb_elem)=-1.0!+2.0*alpha_dummy
       Fp(Dof*nb_elem,1)=0.0
       
       alpha_dummy=-alpha
   
       Fv(1,1)=0.0
       Fv(1,DoF*nb_elem)=0.0

       Fv(Dof*nb_elem,Dof*nb_elem)=0.0
       Fv(Dof*nb_elem,1)=0.0

    end if
  end subroutine init_FpFv


  ! Initializes the diagonal parameters matrices
  subroutine init_DpDv(Dp,Dv,DoF,nb_elem,velocity,density)
    real,dimension(:,:),intent(out) :: Dp
    real,dimension(:,:),intent(out) :: Dv
    integer            ,intent(in)  :: DoF
    integer            ,intent(in)  :: nb_elem
    real,dimension(:)  ,intent(in)  :: velocity
    real,dimension(:)  ,intent(in)  :: density
    integer :: i,j

    Dp=0.0
    Dv=0.0

    do i=1,nb_elem
       do j=1,DoF
          Dp(DoF*(i-1)+j,DoF*(i-1)+j)=velocity(i)**2*density(i)
          Dv(DoF*(i-1)+j,DoF*(i-1)+j)=1.0/density(i)
       end do
    end do
  end subroutine init_DpDv


  ! Initializes the time step length dt
  subroutine init_dt(Ap,Av,dt,nb_elem,time_scheme,velocity)
    real,dimension(:,:),intent(in)  :: Ap
    real,dimension(:,:),intent(in)  :: Av
    real               ,intent(out) :: dt
    integer            ,intent(in)  :: nb_elem
    character(len=*)   ,intent(in)  :: time_scheme
    real,dimension(:)  ,intent(in)  :: velocity
    type(sparse_matrix)             :: sparse_A
    real,dimension(size(Ap,1),size(Ap,2)) :: A
    real                            :: max_value
    integer                         :: i
    real,parameter                  :: alpha=1.80
    real                            :: cfl

    dt=0.0
    
    A=matmul(Ap,Av)
 
    call Full2Sparse(A,sparse_A)
    call power_method_sparse(sparse_A,max_value)
    dt=alpha/(sqrt(abs(max_value)))
    
    call free_sparse_matrix(sparse_A)
    A=matmul(Av,Ap)
    call Full2Sparse(matmul(Av,Ap),sparse_A)
    call power_method_sparse(sparse_A,max_value)
    dt=min(dt,alpha/(sqrt(abs(max_value))))
    call free_sparse_matrix(sparse_A)
    dt=dt*1.0/(maxval(velocity))  !CFL for LF

    if (time_scheme.eq.'LF4') then
       dt=1.0*dt
    else if (time_scheme.eq.'RK4') then
       dt=1.4*dt
    else if (time_scheme.eq.'AB3') then
       dt=dt/2.8
    end if    
  end subroutine init_dt
  

  ! Initializes the main operators (if there is no ABC)
  subroutine init_adjoint_operator(problem)
    type(adjoint_problem),intent(inout)        :: problem
    real,dimension(problem%nb_elem*problem%DoF,                                 &
                   problem%nb_elem*problem%DoF) :: Ap_full,Av_full,B,App_full,  &
                                                   m_glob,minv_glob,s_glob,     &
                                                   Fp,Fv,Dp,Dv
    real,dimension(problem%DoF,problem%DoF)     :: m_loc,minv_loc,s_loc
    integer                                     :: i,j
    integer                                     :: DoF,nb_elem

    type(sparse_matrix)                         :: sparse_dummy
    
    DoF=problem%DoF           ! For sake of readibility
    nb_elem=problem%nb_elem   ! For sake of readibility
    
    call init_quadrature
    call init_m_loc(m_loc,DoF,problem%bernstein)
    call init_m_glob(m_glob,m_loc,nb_elem)
    call init_minv_loc(minv_loc,m_loc)
    call init_minv_glob(minv_glob,minv_loc,nb_elem)
    call init_s_loc(s_loc,DoF,problem%bernstein)
    call init_s_glob(s_glob,s_loc,nb_elem)
    call init_FpFv(Fp,Fv,DoF,nb_elem,problem%boundaries,problem%alpha,          &
                   problem%receiver_loc)
    call init_DpDv(Dp,Dv,DoF,nb_elem,problem%velocity,problem%density)
    
    Ap_full=0.0
    Av_full=0.0
    call Full2Sparse(m_glob,problem%M)
    call Full2Sparse(minv_glob,problem%Minv)

    if (problem%bernstein) then

       Ap_full=(s_glob+matmul(minv_glob,Fp))
       Av_full=(s_glob+matmul(minv_glob,Fv))
       Ap_full=matmul(Dp,Ap_full)
       Av_full=matmul(Dv,Av_full)

       call init_dt((1/problem%dx)*Ap_full,(1/problem%dx)*Av_full,problem%dt,   &
                    problem%nb_elem,problem%time_scheme,problem%velocity)

       if (problem%time_scheme.eq.'LF4') then
          B=matmul(Ap_full,Av_full)
          Av_full=(1.0/24)*(problem%dt/problem%dx)**2*matmul(Av_full,B)+Av_full
          Ap_full=(1.0/24)*(problem%dt/problem%dx)**2*matmul(B,Ap_full)+Ap_full
       end if

       call Full2Sparse(Ap_full,problem%Ap)
       call Full2Sparse(Av_full,problem%Av)

       problem%Ap%Values=(problem%dt/problem%dx)*problem%Ap%Values
       problem%Av%Values=(problem%dt/problem%dx)*problem%Av%Values

       call Full2Sparse(matmul(Dp,minv_glob),problem%Minv_p)
       call Full2Sparse(matmul(Dv,minv_glob),problem%Minv_v)

       problem%Minv_p%Values=(problem%dt/problem%dx)*problem%Minv_p%Values
       problem%Minv_v%Values=(problem%dt/problem%dx)*problem%Minv_v%Values       
    else
       
       Ap_full=s_glob+Fp
       Av_full=s_glob+Fv   

       Ap_full=matmul(minv_glob,Ap_full)
       Ap_full=matmul(Dp,Ap_full)
       Av_full=matmul(minv_glob,Av_full)
       Av_full=matmul(Dv,Av_full)

       call init_dt((1/problem%dx)*Ap_full,(1/problem%dx)*Av_full,problem%dt,   &
                     problem%nb_elem,problem%time_scheme,problem%velocity)

       if (problem%time_scheme.eq.'LF4') then
          B=matmul(Ap_full,Av_full)
          Av_full=(1.0/24)*(problem%dt/problem%dx)**2*matmul(Av_full,B)+Av_full
          Ap_full=(1.0/24)*(problem%dt/problem%dx)**2*matmul(B,Ap_full)+Ap_full
       end if

       call Full2Sparse(Ap_full,problem%Ap)
       call Full2Sparse(Av_full,problem%Av)

       problem%Ap%Values=(problem%dt/problem%dx)*problem%Ap%Values
       problem%Av%Values=(problem%dt/problem%dx)*problem%Av%Values
       
       call Full2Sparse(matmul(Dp,minv_glob),problem%Minv_p)
       call Full2Sparse(matmul(Dv,minv_glob),problem%Minv_v)

       problem%Minv_p%Values=(problem%dt/problem%dx)*problem%Minv_p%Values
       problem%Minv_v%Values=(problem%dt/problem%dx)*problem%Minv_v%Values
    end if

    App_full=0.0
    if (problem%time_scheme.eq.'LF') then
       do i=1,size(App_full,1)
          App_full(i,i)=-1.0
       end do
    end if
    call Full2Sparse(App_full,problem%App)


    if (problem%strategy.eq.'DTA') then
       call init_sparse_matrix(sparse_dummy,problem%Ap%NNN,problem%Ap%nb_ligne, &
            problem%Ap%IA,problem%Ap%JA,problem%Ap%Values)
       call free_sparse_matrix(problem%Ap)
       call transpose_sparse(problem%Av,problem%Ap)
       call free_sparse_matrix(problem%Av)
       call transpose_sparse(sparse_dummy,problem%Av)
       call free_sparse_matrix(sparse_dummy)

       if (problem%scalar_product.eq.'M') then
          problem%App=sparse_matmul(problem%App,problem%M)
          problem%Ap=sparse_matmul(problem%Ap,problem%M)
          problem%Av=sparse_matmul(problem%Av,problem%M)

          problem%App=sparse_matmul(problem%Minv,problem%App)
          problem%Ap=sparse_matmul(problem%Minv,problem%Ap)
          problem%Av=sparse_matmul(problem%Minv,problem%Av)
       end if

    end if
  end subroutine init_adjoint_operator


  ! Initializes the main operators with ABC
  subroutine init_adjoint_operator_abc(problem)
    type(adjoint_problem),intent(inout)        :: problem
    real,dimension(problem%nb_elem*problem%DoF,                                 &
                   problem%nb_elem*problem%DoF) :: Ap_full,Av_full,B,App_full,  &
                                                   m_glob,minv_glob,s_glob,     &
                                                   Fp,Fv,Dp,Dv,dummy,           &
                                                   minv_glob_abc,mabs
    real,dimension(problem%DoF,problem%DoF)     :: m_loc,minv_loc,s_loc
    integer                                     :: i,j
    integer                                     :: DoF,nb_elem
    type(sparse_matrix)                         :: sparse_dummy
    
    DoF=problem%DoF           ! For sake of lisibility
    nb_elem=problem%nb_elem
    
    call init_quadrature
    call init_m_loc(m_loc,DoF,problem%bernstein)
    call init_m_glob(m_glob,m_loc,nb_elem)
    call init_minv_loc(minv_loc,m_loc)
    call init_minv_glob(minv_glob,minv_loc,nb_elem)
    call init_s_loc(s_loc,DoF,problem%bernstein)
    call init_s_glob(s_glob,s_loc,nb_elem)
    call init_FpFv(Fp,Fv,DoF,nb_elem,problem%boundaries,problem%alpha,          &
                   problem%receiver_loc)
    call init_DpDv(Dp,Dv,DoF,nb_elem,problem%velocity,problem%density)

    Ap_full=0.0
    Av_full=0.0
    App_full=0.0
    call Full2Sparse(m_glob,problem%M)
    call Full2Sparse(minv_glob,problem%Minv)
    
    if (problem%bernstein) then
       
       Av_full=(s_glob+matmul(minv_glob,Fv))
       Av_full=matmul(Dv,Av_full)



       call init_dt((1/problem%dx)*Av_full,(1/problem%dx)*Av_full,problem%dt,   &
                    problem%nb_elem,problem%time_scheme,problem%velocity)

       call init_minv_glob_abc(minv_glob_abc,m_loc,nb_elem,DoF,                 &
                               problem%velocity(1),problem%velocity(nb_elem),   &
                               problem%dx,problem%dt,problem%time_scheme)
       call init_mabs(mabs,nb_elem,DoF,                                         &
                      problem%velocity(1),problem%velocity(nb_elem),            &
                      problem%dx,problem%dt,problem%time_scheme)
       
       Ap_full=matmul(m_glob,s_glob)+Fp
       Ap_full=matmul(Dp,Ap_full)
       Ap_full=matmul(minv_glob_abc,Ap_full)


       if (problem%time_scheme.eq.'LF') then
          App_full=mabs-m_glob
       else
          App_full=mabs
       end if
       App_full=matmul(minv_glob_abc,App_full)

       if (problem%time_scheme.eq.'LF4') then
          B=matmul(Ap_full,Av_full)
          Av_full=(1.0/24)*(problem%dt/problem%dx)**2*matmul(Av_full,B)+Av_full
          Ap_full=(1.0/24)*(problem%dt/problem%dx)**2*matmul(B,Ap_full)+Ap_full
       end if
       

          call Full2Sparse(Ap_full,problem%Ap)
          call Full2Sparse(Av_full,problem%Av)
          call Full2Sparse(App_full,problem%App)
       
       problem%Ap%Values=(problem%dt/problem%dx)*problem%Ap%Values
       problem%Av%Values=(problem%dt/problem%dx)*problem%Av%Values

       call Full2Sparse(matmul(Dp,minv_glob_abc),problem%Minv_p)
       call Full2Sparse(matmul(Dv,minv_glob),problem%Minv_v)
       
       problem%Minv_p%Values=(problem%dt/problem%dx)*problem%Minv_p%Values
       problem%Minv_v%Values=(problem%dt/problem%dx)*problem%Minv_v%Values
    else
       
       Av_full=s_glob+Fv   


       Av_full=matmul(minv_glob,Av_full)
       Av_full=matmul(Dv,Av_full)

       call init_dt((1/problem%dx)*Av_full,(1/problem%dx)*Av_full,problem%dt,   &
                     problem%nb_elem,problem%time_scheme,problem%velocity)


       call init_minv_glob_abc(minv_glob_abc,m_loc,nb_elem,DoF,                 &
                               problem%velocity(1),problem%velocity(nb_elem),   &
                               problem%dx,problem%dt,problem%time_scheme)
       call init_mabs(mabs,nb_elem,DoF,problem%velocity(1),                     &
                      problem%velocity(nb_elem),problem%dx,problem%dt,          &
                      problem%time_scheme)
       
       Ap_full=s_glob+Fp
       
       Ap_full=matmul(minv_glob_abc,Ap_full)
       Ap_full=matmul(Dp,Ap_full)
       
       if (problem%time_scheme.eq.'LF') then
          App_full=mabs-m_glob
       else
          App_full=mabs
       end if
       App_full=matmul(minv_glob_abc,App_full)

       if (problem%time_scheme.eq.'LF4') then
          B=matmul(Ap_full,Av_full)
          Av_full=(1.0/24)*(problem%dt/problem%dx)**2*matmul(Av_full,B)+Av_full
          Ap_full=(1.0/24)*(problem%dt/problem%dx)**2*matmul(B,Ap_full)+Ap_full
       end if

          call Full2Sparse(Ap_full,problem%Ap)
          call Full2Sparse(Av_full,problem%Av)
          call Full2Sparse(App_full,problem%App)

       problem%Ap%Values=(problem%dt/problem%dx)*problem%Ap%Values
       problem%Av%Values=(problem%dt/problem%dx)*problem%Av%Values

       call Full2Sparse(matmul(Dp,minv_glob_abc),problem%Minv_p)
       call Full2Sparse(matmul(Dv,minv_glob),problem%Minv_v)

       problem%Minv_p%Values=(problem%dt/problem%dx)*problem%Minv_p%Values
       problem%Minv_v%Values=(problem%dt/problem%dx)*problem%Minv_v%Values
    end if

    if (problem%strategy.eq.'DTA') then
       call init_sparse_matrix(sparse_dummy,problem%Ap%NNN,problem%Ap%nb_ligne, &
            problem%Ap%IA,problem%Ap%JA,problem%Ap%Values)
       call free_sparse_matrix(problem%Ap)
       call transpose_sparse(problem%Av,problem%Ap)
       call free_sparse_matrix(problem%Av)
       call transpose_sparse(sparse_dummy,problem%Av)
       call free_sparse_matrix(sparse_dummy)
       call init_sparse_matrix(sparse_dummy,problem%App%NNN,                    &
                               problem%App%nb_ligne,problem%App%IA,             &
                               problem%App%JA,problem%App%Values)
       call free_sparse_matrix(problem%App)
       call transpose_sparse(sparse_dummy,problem%App)
       call free_sparse_matrix(sparse_dummy)

       if (problem%scalar_product.eq.'M') then
          problem%App=sparse_matmul(problem%App,problem%M)
          problem%Ap=sparse_matmul(problem%Ap,problem%M)
          problem%Av=sparse_matmul(problem%Av,problem%M)

          problem%App=sparse_matmul(problem%Minv,problem%App)
          problem%Ap=sparse_matmul(problem%Minv,problem%Ap)
          problem%Av=sparse_matmul(problem%Minv,problem%Av)
       end if
    end if
    
  end subroutine init_adjoint_operator_abc
  
  
  !**************** PROBLEM RESOLUTION ******************************************
  ! Makes one backward time step
  subroutine one_backward_step(problem,t)
    type(adjoint_problem),intent(inout) :: problem
    real,                  intent(in)    :: t

    real,dimension(problem%DoF*problem%nb_elem) :: RHSu0,RHSp0
    real,dimension(problem%DoF*problem%nb_elem) :: RHSu1,RHSp1
    real,dimension(problem%DoF*problem%nb_elem) :: RHSu2,RHSp2
    real,dimension(problem%DoF*problem%nb_elem) :: U1,U2,U3,U4
    real,dimension(problem%DoF*problem%nb_elem) :: P1,P2,P3,P4
    real                                        :: gd,gg
    real                                        :: f0
    integer                                     :: i
    real,dimension(problem%DoF*problem%nb_elem) :: FP_0,FP_half,FP_1
    real,dimension(problem%DoF*problem%nb_elem) :: FU_0,FU_half,FU_1

    real,dimension(problem%nb_elem*problem%DoF) :: zero

    zero=0.0

    if ((problem%time_scheme.eq.'LF').or.(problem%time_scheme.eq.'LF4')) then

       call eval_RHS(FP_0,FU_0,problem,t-problem%dt)
       call eval_RHS(FP_half,FU_half,problem,t-0.5*problem%dt)

       call LF_forward(problem%P,problem%U,problem%Ap,problem%Av,problem%App,   &
            FP_0,FU_half)

    else if (problem%time_scheme.eq.'RK4') then

       call eval_RHS(FP_0,FU_0,problem,t-problem%dt)

       call eval_RHS(FP_half,FU_half,problem,t-0.5*problem%dt)
       problem%RHSu_half=FU_half
       problem%RHSp_half=FP_half
       call eval_RHS(FP_1,FU_1,problem,t)
       problem%RHSu=FU_1
       problem%RHSp=FP_1

       if (problem%strategy.eq.'DTA') then
          call RK4_forward(problem%P,problem%U,problem%Ap,problem%Av,problem%App,  &
               zero,zero,zero,zero,zero,zero,FP_1,FU_1)
       else
          call RK4_forward(problem%P,problem%U,problem%Ap,problem%Av,problem%App,  &
               FP_0,FP_half,FP_1,FU_0,FU_half,FU_1)
       end if

    else if (problem%time_scheme.eq.'AB3') then

       RHSp0=0.0
       RHSu0=0.0
       RHSp1=0.0
       RHSu1=0.0
       RHSp2=0.0
       RHSu2=0.0

       if (problem%strategy.eq.'DTA') then
          call eval_RHS(problem%RHSp,problem%RHSu,problem,t-problem%dt)
          call AB3_forward(problem%P,problem%U,problem%Ap,problem%Av,problem%App,  &
               zero,zero,                                              &
               zero,zero,                                              &
               zero,zero,                                              &
               problem%Pk1,problem%Pk2,problem%Uk1,problem%Uk2,        &
               problem%RHSp,problem%RHSu)
       else
          call eval_RHS(problem%RHSp,problem%RHSu,problem,t-problem%dt)
          call eval_RHS(RHSp0,RHSu0,problem,t-problem%dt)
          call eval_RHS(RHSp2,RHSu1,problem,t-2*problem%dt)
          call eval_RHS(RHSp2,RHSu2,problem,t-3*problem%dt)
          call AB3_forward(problem%P,problem%U,problem%Ap,problem%Av,problem%App,  &
               RHSp0,RHSu0,                                              &
               RHSp1,RHSu1,                                              &
               RHSp2,RHSu2,                                              &
               problem%Pk1,problem%Pk2,problem%Uk1,problem%Uk2)
       end if

    end if
  end subroutine one_backward_step


  ! Makes all the backward steps
  subroutine all_backward_step(problem)
    type(adjoint_problem),intent(inout) :: problem
    integer                              :: i
    real                                 :: t

    do i=1,problem%n_time_step
       t=i*problem%dt
       call one_backward_step(problem,t)
    end do    
  end subroutine all_backward_step


  ! Frees the acoustic adjoint problem variables
  subroutine free_adjoint_problem(problem)
    type(adjoint_problem),intent(inout) :: problem
    deallocate(problem%U)
    deallocate(problem%P)
    deallocate(problem%density)
    deallocate(problem%velocity)
    deallocate(problem%data_P)
    deallocate(problem%data_U)
    deallocate(problem%P_received)
    deallocate(problem%U_received)
    deallocate(problem%RHSu)
    deallocate(problem%RHSp)
    if (problem%time_scheme.eq.'RK4') then
    deallocate(problem%RHSu_half)
    deallocate(problem%RHSp_half)
    end if
    call free_sparse_matrix(problem%Ap)
    call free_sparse_matrix(problem%Av)
    call free_sparse_matrix(problem%Minv_p)
    call free_sparse_matrix(problem%Minv_v)
    call free_sparse_matrix(problem%M)
    call free_sparse_matrix(problem%Minv)
    call free_sparse_matrix(problem%App)
  end subroutine free_adjoint_problem



  ! Evaluates the backward source (here the source are the receiver)
  function eval_backward_source(P_received,data_P,t,final_time)
    real,dimension(:,:),intent(in) :: P_received
    real,dimension(:,:),intent(in) :: data_P
    real               ,intent(in) :: t
    real               ,intent(in) :: final_time
    real                           :: eval_P_received
    real                           :: eval_data_P
    real                           :: eval_backward_source
    integer                        :: i
    real                           :: t2

    t2=t
    i=1
    do while ((P_received(i,1).le.t2).and.(i.ne.size(data_P,1)))
       i=i+1
    end do
    i=i-1

    if (i.eq.0) then
       print*,'test',i,t,t2,P_received(1,1)
    end if

    eval_P_received=(t2-P_received(i,1))*P_received(i+1,2)       &
                        +(P_received(i+1,1)-t2)*P_received(i,2)
    eval_P_received=eval_P_received/(P_received(i+1,1)           &
         -P_received(i,1))

    i=2
    do while ((data_P(i,1).le.t2).and.(i.ne.size(data_P,1)))
       i=i+1
    end do
    i=i-1

    if (i.eq.0) then
       print*,'test',i,t,t2
    end if

    eval_data_P=(t2-data_P(i,1))*data_P(i+1,2)       &
                        +(data_P(i+1,1)-t2)*data_P(i,2)
    eval_data_P=eval_data_P/(data_P(i+1,1)           &
         -data_P(i,1))


    eval_backward_source=eval_data_P-eval_P_received
  end function eval_backward_source


  ! Evaluates the RHS for the Backward problem at time=t
  subroutine eval_RHS(RHSp,RHSu,problem,t)
    real,dimension(:)   ,intent(inout) :: RHSp
    real,dimension(:)   ,intent(inout) :: RHSu
    type(adjoint_problem),intent(in)  :: problem
    real                  ,intent(in)  :: t
    
    real    :: gg,gd,f0
    integer :: last_node
    
    RHSp=0.0
    RHSu=0.0
    f0=2.5

    last_node=problem%DoF*problem%nb_elem

    if (t.le.0.0) then
       RHSp=0.0
       RHSu=0.0
    else

       if (problem%boundaries.eq.'dirichlet') then
          gg=min(0.5*t,0.5)
          gd=-gg
          RHSu(1)=gg*(-1.0+0.0*problem%alpha)
          RHSu(last_node)=gd*(1.0-0.0*problem%alpha)
       end if

       if (problem%signal.eq.'flat') then
         ! gg=(2*t-1/f0)*exp(-(2.0*PI*(2*t-1/f0)*f0)**2.0)*5/0.341238111
         ! RHSu((problem%source_loc-1)*problem%DoF+1)=gg*(-1.0-0.0*problem%alpha)
          
          gg=eval_backward_source(problem%P_received,problem%data_P,t,problem%final_time)
          !  RHSu(problem%receiver_loc*problem%DoF+1)=gg*(-1.0-2.0*problem%alpha) 
          RHSp((problem%receiver_loc-1)*problem%DoF+1)=gg*(1.0-0.0*problem%alpha) 
       end if
    end if

    if (problem%strategy.eq.'ATD') then
       RHSp=sparse_matmul(problem%Minv_p,RHSp)
       RHSu=sparse_matmul(problem%Minv_v,RHSu)
    end if

  end subroutine eval_RHS
end module m_adjoint
