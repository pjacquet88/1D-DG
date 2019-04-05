module m_acoustic
  use m_polynom
  use m_matrix
  use m_powermethod
  use m_file_function
  use m_time_scheme
  implicit none


  ! Acoustic problem type
  type acoustic_problem
     integer                            :: DoF          ! nb of degree of freedom
     integer                            :: nb_elem      ! nb of elements
     character(len=20)                  :: time_scheme  ! the used time scheme
     logical                            :: bernstein    ! bool for bernstein elmnts
     real(mp),dimension(:) ,allocatable :: density      ! density model
     real(mp),dimension(:) ,allocatable :: velocity     ! velocity model
     real(mp)                           :: dx           ! element length
     real(mp)                           :: dt           ! step time length
     real(mp)                           :: total_length ! domain length
     real(mp)                           :: final_time   ! final time
     real(mp)                           :: alpha        ! penalisation value
     real(mp),dimension(:),allocatable  :: U            ! velocity unknowns vector
     real(mp),dimension(:),allocatable  :: P            ! pressure unknowns vector
     real(mp),dimension(:),allocatable  :: Pk1,Pk2      ! AB3 vectors
     real(mp),dimension(:),allocatable  :: Uk1,Uk2      ! AB3 vectors
     character(len=20)                  :: boundaries   ! string for BC
     character(len=20)                  :: signal       ! string to initialize U,P
     character(len=20)                  :: flux         ! Flux for DG
     real(mp),dimension(:),allocatable  :: RHSu
     real(mp),dimension(:),allocatable  :: RHSp
     real(mp),dimension(:),allocatable  :: RHSu_half
     real(mp),dimension(:),allocatable  :: RHSp_half

     !----------- MATRICES ---------------------------
     type(sparse_matrix)             :: Ap,App
     type(sparse_matrix)             :: Av
     type(sparse_matrix)             :: Minv_p,Minv_v
     type(sparse_matrix)             :: M,Minv
     type(sparse_matrix)             :: S
     !----------- SOURCE & RECEIVER ------------------
     integer                         :: source_loc    ! beginning of an element
     integer                         :: receiver_loc  ! beginning of an element
     integer                         :: n_time_step   ! nb of time step
  end type acoustic_problem

  real(mp),dimension(15) :: xx,weight
  real(mp),parameter     :: PI=acos(-1.0_mp)

  public  :: init_acoustic_problem,init_acoustic_operator,all_forward_step,     &
             one_forward_step,                                                  &
             free_acoustic_problem,                                             &
             error_periodic

  private :: xx,weight,initial_perturbation,PI,                                  &
             init_quadrature,Int_bxb,Int_lxl,init_m_loc,init_s_loc,             &
             init_UP,init_m_glob,init_minv_glob,init_minv_glob_abc,init_mabs,   &
             init_s_glob,init_FpFv,init_DpDv,init_dt
  
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
       xx(i)=1.0_mp-xx(16-i)
       weight(i)=weight(16-i)
    end do
  end subroutine init_quadrature

  ! Evaluates the integral of the product of 2 bernstein polyniomials
  function Int_bxb(bpol1,bpol2)
    real(mp),dimension(:),intent(in) :: bpol1
    real(mp),dimension(:),intent(in) :: bpol2
    real(mp)                         :: Int_bxb
    integer                           :: i

    Int_bxb=0.0_mp

    do i=1,15
       Int_bxb=Int_bxb                  &
            +eval_polynom_b(bpol1,xx(i))&
            *eval_polynom_b(bpol2,xx(i))&
            *weight(i)
    end do

  end function Int_bxb

  ! Evaluates the integral of the product of 2 lagrange polynomials  
  function Int_lxl(lpol1,lpol2)
    real(mp),dimension(:),intent(in) :: lpol1
    real(mp),dimension(:),intent(in) :: lpol2
    real(mp)                         :: Int_lxl
    integer                          :: i

    Int_lxl=0.0_mp

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
    real(mp),intent(in) :: x
    real(mp),intent(in) :: dt
    character(len=*),intent(in) :: signal
    real(mp)                    :: initial_perturbation 
    real(mp)                    :: xt

    xt=x+dt
    
    if (signal.eq.'sinus') then
       initial_perturbation=sin(4*PI*xt)
       ! initial_perturbation=sin(xt)
    else if (signal.eq.'boxcar') then
       if ((x.lt.0.55_mp).and.(x.gt.0.45_mp)) then
          initial_perturbation=0.5_mp
       else
          initial_perturbation=0.0_mp
       end if
    else if (signal.eq.'x') then
       initial_perturbation=xt**5
        else if (signal.eq.'flat') then
       initial_perturbation=0.0_mp
    else
       initial_perturbation=20.0_mp*(xt-0.5_mp)                           &
                            *exp(-(2.0_mp*PI*(xt-0.5_mp)*2.0_mp)**2.0_mp)
    end if    
  end function initial_perturbation


  ! Initializes the acoustic problem
  subroutine init_acoustic_problem(problem,nb_elem,DoF,time_scheme,velocity,    &
                                   density,total_length,final_time,alpha,       &
                                   bernstein,flux,signal, boundaries,source_loc,&
                                   receiver_loc)
    
    type(acoustic_problem),intent(out) :: problem
    integer               ,intent(in)  :: nb_elem
    integer               ,intent(in)  :: DoF
    real(mp),dimension(:) ,intent(in)  :: velocity
    real(mp),dimension(:) ,intent(in)  :: density
    character(len=*)      ,intent(in)  :: time_scheme
    real(mp)              ,intent(in)  :: total_length
    real(mp)              ,intent(in)  :: final_time
    real(mp)              ,intent(in)  :: alpha
    logical               ,intent(in)  :: bernstein
    character(len=*)      ,intent(in)  :: flux
    character(len=*)      ,intent(in)  :: signal
    character(len=*)      ,intent(in)  :: boundaries
    integer               ,intent(in)  :: source_loc
    integer               ,intent(in)  :: receiver_loc
    
    integer  :: i,j
    integer  :: dv,dp,size_velocity,size_density
    real(mp) :: gd,dd

    problem%DoF=DoF
    problem%time_scheme=time_scheme
    problem%nb_elem=nb_elem
    problem%total_length=total_length
    problem%dx=total_length/(nb_elem)
    problem%final_time=final_time
    problem%alpha=alpha
    problem%bernstein=bernstein
    problem%boundaries=boundaries
    problem%flux=flux
    problem%signal=signal
    problem%receiver_loc=receiver_loc
    problem%source_loc=source_loc
    
    size_velocity=size(velocity)
    size_density=size(density)
    
    allocate(problem%U(DoF*nb_elem),problem%P(DoF*nb_elem))
    problem%U=0.0_mp
    problem%P=0.0_mp
    
    if (problem%time_scheme.eq.'AB3') then
       allocate(problem%Pk1(DoF*nb_elem),problem%Uk1(DoF*nb_elem))
       allocate(problem%Pk2(DoF*nb_elem),problem%Uk2(DoF*nb_elem))
       problem%Pk1=0.0_mp
       problem%Pk2=0.0_mp
       problem%Uk1=0.0_mp
       problem%Uk2=0.0_mp
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


    !---------------init matrices------------------
    if (problem%boundaries.eq.'ABC') then
       call init_acoustic_operator_abc(problem)
    else
       call init_acoustic_operator(problem)
    end if

    problem%n_time_step=int(problem%final_time/problem%dt)+1
    problem%final_time=problem%dt*problem%n_time_step
    ! here I change the final time to fit frames of the forward and the backward

    allocate(problem%RHSu(nb_elem*DoF))
    allocate(problem%RHSp(nb_elem*DoF))
    if (time_scheme.eq.'RK4') then
       allocate(problem%RHSu_half(nb_elem*DoF))
       allocate(problem%RHSp_half(nb_elem*DoF))
     end if

     call init_UP(problem)

  end subroutine init_acoustic_problem

  ! Initialize P and U at time=0  
  subroutine init_UP(problem)
    type(acoustic_problem),intent(inout) :: problem
    integer                              :: i,j,DoF
    real(mp)                             :: x,ddx,dephasing
    
    DoF=problem%DoF
    ddx=problem%dx/(DoF-1)
    x=0.0_mp
    
    if ((problem%time_scheme.eq.'LF').or.(problem%time_scheme.eq.'LF4')) then
       dephasing=-problem%dt/2.0_mp
    else
       dephasing=0.0_mp
    end if

    do i=1,problem%nb_elem
      do j=1,DoF
          x=(i-1)*problem%dx+(j-1)*ddx
          problem%U((i-1)*problem%DoF+j)=initial_perturbation(x,0.0_mp,problem%signal)
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
    real(mp),dimension(:,:),intent(out) :: m_loc
    integer                ,intent(in)  :: DoF
    logical                ,intent(in)  :: bernstein
    integer                             :: i,j

    m_loc=0.0_mp

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
    real(mp),dimension(:,:),intent(out) :: minv_loc
    real(mp),dimension(:,:),intent(in)  :: m_loc
    minv_loc=0.0_mp
    minv_loc=LU_inv(m_loc)
  end subroutine init_minv_loc


  ! Initialize the global mass matrix (on all the elements)
  subroutine init_m_glob(m_glob,m_loc,nb_elem)
    real(mp),dimension(:,:),intent(out) :: m_glob
    real(mp),dimension(:,:),intent(in)  :: m_loc
    integer                ,intent(in)  :: nb_elem

    integer :: i,DoF
    m_glob=0.0_mp
    DoF=size(m_loc,1)
        
    do i=1,nb_elem
       m_glob(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=m_loc
    end do
  end subroutine init_m_glob
  

  ! Initialize the global inverse mass matrix (on all the element)
  subroutine init_minv_glob(minv_glob,minv_loc,nb_elem)
    real(mp),dimension(:,:),intent(out) :: minv_glob
    real(mp),dimension(:,:),intent(in)  :: minv_loc
    integer                ,intent(in)  :: nb_elem
    integer                             :: i,DoF

    minv_glob=0.0_mp

    DoF=size(minv_loc,1)

    do i=1,nb_elem
       minv_glob(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=minv_loc
    end do
  end subroutine init_minv_glob


  ! initialize the global inverse "abc mass matrix" 
  subroutine init_minv_glob_abc(minv_glob_abc,m_loc,nb_elem,DoF,velocity_beg,   &
                                velocity_end,dx,dt,time_scheme)
    real(mp),dimension(:,:),intent(out) :: minv_glob_abc
    real(mp),dimension(:,:),intent(in)  :: m_loc
    integer                ,intent(in)  :: nb_elem
    integer                ,intent(in)  :: DoF
    real(mp)               ,intent(in)  :: velocity_beg
    real(mp)               ,intent(in)  :: velocity_end
    real(mp)               ,intent(in)  :: dx
    real(mp)               ,intent(in)  :: dt
    character(len=*)       ,intent(in)  :: time_scheme

    integer                     :: i
    real(mp),dimension(DoF,DoF) :: m_loc_abc1,m_loc_abc2
    real(mp),dimension(DoF,DoF) :: minv_loc,minv_loc_abc1,minv_loc_abc2

    minv_glob_abc=0.0_mp

    m_loc_abc1=m_loc
    m_loc_abc2=m_loc

    if (time_scheme.eq.'LF') then
       m_loc_abc2(DoF,DoF)=m_loc_abc2(DoF,DoF)+0.5_mp*velocity_end*dt/dx
       m_loc_abc1(1,1)=m_loc_abc1(1,1)+0.5_mp*velocity_beg*dt/dx
    end if

    minv_loc=LU_inv(m_loc)
    minv_loc_abc1=LU_inv(m_loc_abc1)
    minv_loc_abc2=LU_inv(m_loc_abc2)

    minv_glob_abc=0.0_mp
    minv_glob_abc(1:DoF,1:DoF)=minv_loc_abc1
    do i=2,nb_elem-1
       minv_glob_abc(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=minv_loc
    end do
    minv_glob_abc(DoF*(nb_elem-1)+1:DoF*nb_elem,DoF*(nb_elem-1)+1:DoF*nb_elem)  &
         =minv_loc_abc2
  end subroutine init_minv_glob_abc


  ! Initializes the global Abs matrix
  subroutine init_mabs(mabs,nb_elem,DoF,velocity_beg,velocity_end,dx,dt,time_scheme)
    real(mp),dimension(:,:),intent(out) :: mabs
    integer                ,intent(in)  :: nb_elem
    integer                ,intent(in)  :: DoF
    real(mp)               ,intent(in)  :: velocity_beg
    real(mp)               ,intent(in)  :: velocity_end
    real(mp)               ,intent(in)  :: dx
    real(mp)               ,intent(in)  :: dt
    character(len=*)       ,intent(in)  :: time_scheme

    integer :: i

    mabs=0.0_mp

    if (time_scheme.eq.'LF') then
       mabs(DoF*nb_elem,DoF*nb_elem)=0.5_mp*velocity_end*dt/dx
       mabs(1,1)=0.5_mp*velocity_beg*dt/dx
    else
       mabs(DoF*nb_elem,DoF*nb_elem)=velocity_end*dt/dx
       mabs(1,1)=velocity_beg*dt/dx
    end if
  end subroutine init_mabs


  ! Initailizes the local stiffness matrix (on one element)
  subroutine init_s_loc(s_loc,DoF,bernstein)
    real(mp),dimension(:,:),intent(out) :: s_loc
    integer                ,intent(in)  :: DoF
    logical                ,intent(in)  :: bernstein
    integer                             :: i,j
    real(mp),dimension(DoF)             :: bpol,dbpol
    real(mp),dimension(DoF)             :: dlpol

    s_loc=0.0_mp
    
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
    real(mp),dimension(:,:),intent(out) :: s_glob
    real(mp),dimension(:,:),intent(in)  :: s_loc
    integer                ,intent(in)  :: nb_elem
    integer :: i,DoF

    s_glob=0.0_mp
    DoF=size(s_loc,1)


    do i=1,nb_elem
       s_glob(DoF*(i-1)+1:Dof*i,DoF*(i-1)+1:Dof*i)=s_loc
    end do
  end subroutine init_s_glob


  ! Initalize the Fp and Fv flux matrices (with penalization)
  subroutine init_FpFv(Fp,Fv,DoF,nb_elem,boundaries,flux,alpha,receiver_loc)
    real(mp),dimension(:,:),intent(out) :: Fp
    real(mp),dimension(:,:),intent(out) :: Fv
    integer                ,intent(in)  :: DoF
    integer                ,intent(in)  :: nb_elem
    character(len=*)       ,intent(in)  :: boundaries
    character(len=*)       ,intent(in)  :: flux
    real(mp)               ,intent(in)  :: alpha
    integer                ,intent(in)  :: receiver_loc

    integer  :: i,j
    real(mp) :: alpha_dummy

    Fp=0.0_mp
    Fv=0.0_mp

    select case (flux)
    case('center')

      alpha_dummy=alpha

      Fp(DoF,DoF)=-0.5_mp+alpha_dummy
      Fp(DoF,DoF+1)=0.5_mp-alpha_dummy

      do i=2,nb_elem-1
        Fp(Dof*(i-1)+1,Dof*(i-1)+1)=0.5_mp+alpha_dummy
        Fp(Dof*(i-1)+1,Dof*(i-1))=-0.5_mp-alpha_dummy

        Fp(Dof*i,Dof*i)=-0.5_mp+alpha_dummy
        Fp(Dof*i,Dof*i+1)=0.5_mp-alpha_dummy
      end do
      Fp(Dof*(nb_elem-1)+1,Dof*(nb_elem-1)+1)=0.5_mp+alpha_dummy
      Fp(Dof*(nb_elem-1)+1,Dof*(nb_elem-1))=-0.5_mp-alpha_dummy

      alpha_dummy=-alpha

      Fv(DoF,DoF)=-0.5_mp+alpha_dummy
      Fv(DoF,DoF+1)=0.5_mp-alpha_dummy

      do i=2,nb_elem-1
        Fv(Dof*(i-1)+1,Dof*(i-1)+1)=0.5_mp+alpha_dummy
        Fv(Dof*(i-1)+1,Dof*(i-1))=-0.5_mp-alpha_dummy

        Fv(Dof*i,Dof*i)=-0.5_mp+alpha_dummy
        Fv(Dof*i,Dof*i+1)=0.5_mp-alpha_dummy
      end do
      Fv(Dof*(nb_elem-1)+1,Dof*(nb_elem-1)+1)=0.5_mp+alpha_dummy
      Fv(Dof*(nb_elem-1)+1,Dof*(nb_elem-1))=-0.5_mp-alpha_dummy

      select case (boundaries)

      case ('periodic')

        alpha_dummy=alpha

        Fp(1,1)=0.5_mp+alpha_dummy
        Fp(1,DoF*nb_elem)=-0.5_mp-alpha_dummy

        Fp(Dof*nb_elem,Dof*nb_elem)=-0.5_mp+alpha_dummy
        Fp(Dof*nb_elem,1)=0.5_mp-alpha_dummy

        alpha_dummy=-alpha

        Fv(1,1)=0.5_mp+alpha_dummy
        Fv(1,DoF*nb_elem)=-0.5_mp-alpha_dummy

        Fv(Dof*nb_elem,Dof*nb_elem)=-0.5_mp+alpha_dummy
        Fv(Dof*nb_elem,1)=0.5_mp-alpha_dummy

      case ('dirichlet') 

        alpha_dummy=alpha

        Fp(1,1)=0.0_mp
        Fp(1,DoF*nb_elem)=0.0_mp

        Fp(Dof*nb_elem,Dof*nb_elem)=0.0_mp
        Fp(Dof*nb_elem,1)=0.0_mp


        alpha_dummy=-alpha

        Fv(1,1)=1.0_mp!+2.0_mp*alpha_dummy
        Fv(1,DoF*nb_elem)=0.0_mp

        Fv(Dof*nb_elem,Dof*nb_elem)=-1.0_mp!-2.0_mp*alpha_dummy
        Fv(Dof*nb_elem,1)=0.0_mp


      case ('ABC')

        alpha_dummy=alpha

        Fv(1,1)=0.0_mp!(1.0_mp+2.0_mp*alpha_dummy)
        Fv(1,DoF*nb_elem)=0.0_mp

        Fv(Dof*nb_elem,Dof*nb_elem)=0.0_mp!-1.0_mp-2.0_mp*alpha_dummy
        Fv(Dof*nb_elem,1)=0.0_mp

        alpha_dummy=-alpha

        Fp(1,1)=1.0_mp
        Fp(1,DoF*nb_elem)=0.0_mp

        Fp(Dof*nb_elem,Dof*nb_elem)=-1.0_mp!*2
        Fp(Dof*nb_elem,1)=0.0_mp

      case ('neumann')

        alpha_dummy=alpha

        Fp(1,1)=1.0_mp!-2.0_mp*alpha_dummy
        Fp(1,DoF*nb_elem)=0.0_mp

        Fp(Dof*nb_elem,Dof*nb_elem)=-1.0_mp!+2.0_mp*alpha_dummy
        Fp(Dof*nb_elem,1)=0.0_mp

        alpha_dummy=-alpha

        Fv(1,1)=0.0_mp
        Fv(1,DoF*nb_elem)=0.0_mp

        Fv(Dof*nb_elem,Dof*nb_elem)=0.0_mp
        Fv(Dof*nb_elem,1)=0.0_mp

      end select

    case DEFAULT

      print*,'This flux is not implemented yet'
      STOP

    end select

  end subroutine init_FpFv

! Initializes the diagonal parameter matrices
subroutine init_DpDv(Dp,Dv,DoF,nb_elem,velocity,density)
  real(mp),dimension(:,:),intent(out) :: Dp
  real(mp),dimension(:,:),intent(out) :: Dv
  integer                ,intent(in)  :: DoF
  integer                ,intent(in)  :: nb_elem
  real(mp),dimension(:)  ,intent(in)  :: velocity
  real(mp),dimension(:)  ,intent(in)  :: density
  integer :: i,j

  Dp=0.0_mp
  Dv=0.0_mp

  do i=1,nb_elem
    do j=1,DoF
      Dp(DoF*(i-1)+j,DoF*(i-1)+j)=velocity(i)**2*density(i)
      Dv(DoF*(i-1)+j,DoF*(i-1)+j)=1.0_mp/density(i)
    end do
  end do
end subroutine init_DpDv


  ! Initializes the time step length dt
  subroutine init_dt(Ap,Av,dt,nb_elem,time_scheme,velocity)
    real(mp),dimension(:,:),intent(in)  :: Ap
    real(mp),dimension(:,:),intent(in)  :: Av
    real(mp)               ,intent(out) :: dt
    integer                ,intent(in)  :: nb_elem
    character(len=*)       ,intent(in)  :: time_scheme
    real(mp),dimension(:)  ,intent(in)  :: velocity

    type(sparse_matrix)                       :: sparse_A
    real(mp),dimension(size(Ap,1),size(Ap,2)) :: A
    real(mp)                                  :: max_value
    integer                                   :: i
    real(mp),parameter                        :: alpha=1.80_mp
    real(mp)                                  :: cfl

    dt=0.0_mp

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
    dt=dt*1.0_mp/(maxval(velocity))  !CFL for LF

    if (time_scheme.eq.'LF4') then
       dt=1.0_mp*dt
    else if (time_scheme.eq.'RK4') then
       dt=1.4_mp*dt
    else if (time_scheme.eq.'AB3') then
       dt=dt/2.8_mp
    end if
  end subroutine init_dt

  ! Initializes the main operators (if there is no ABC)
  subroutine init_acoustic_operator(problem)
    type(acoustic_problem),intent(inout)        :: problem
    real(mp),dimension(problem%nb_elem*problem%DoF,                              &
                   problem%nb_elem*problem%DoF) :: Ap_full,Av_full,B,App_full,  &
                                                   m_glob,minv_glob,s_glob,     &
                                                   Fp,Fv,Dp,Dv
    real(mp),dimension(problem%DoF,problem%DoF) :: m_loc,minv_loc,s_loc
    integer                                     :: i,j
    integer                                     :: DoF,nb_elem

    DoF=problem%DoF           ! For sake of readibility
    nb_elem=problem%nb_elem   ! For sake of readibility

    call init_quadrature
    call init_m_loc(m_loc,DoF,problem%bernstein)
    call init_m_glob(m_glob,m_loc,nb_elem)
    call init_minv_loc(minv_loc,m_loc)
    call init_minv_glob(minv_glob,minv_loc,nb_elem)
    call init_s_loc(s_loc,DoF,problem%bernstein)
    call init_s_glob(s_glob,s_loc,nb_elem)
    call init_FpFv(Fp,Fv,DoF,nb_elem,problem%boundaries,problem%flux,           &
                   problem%alpha,problem%receiver_loc)
    call init_DpDv(Dp,Dv,DoF,nb_elem,problem%velocity,problem%density)

    Ap_full=0.0_mp
    Av_full=0.0_mp
    call Full2Sparse(m_glob,problem%M)
    call Full2Sparse(minv_glob,problem%Minv)
    call Full2Sparse(s_glob,problem%S)

    if (problem%bernstein) then

       Ap_full=(s_glob+matmul(minv_glob,Fp))
       Av_full=(s_glob+matmul(minv_glob,Fv))
       Ap_full=matmul(Dp,Ap_full)
       Av_full=matmul(Dv,Av_full)


       call init_dt((1/problem%dx)*Ap_full,(1/problem%dx)*Av_full,problem%dt,   &
                    problem%nb_elem,problem%time_scheme,problem%velocity)

       if (problem%time_scheme.eq.'LF4') then
          B=matmul(Ap_full,Av_full)
          Av_full=(1.0_mp/24)*(problem%dt/problem%dx)**2*matmul(Av_full,B)+Av_full
          Ap_full=(1.0_mp/24)*(problem%dt/problem%dx)**2*matmul(B,Ap_full)+Ap_full
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
          Av_full=(1.0_mp/24)*(problem%dt/problem%dx)**2*matmul(Av_full,B)+Av_full
          Ap_full=(1.0_mp/24)*(problem%dt/problem%dx)**2*matmul(B,Ap_full)+Ap_full
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

    App_full=0.0_mp
    if (problem%time_scheme.eq.'LF') then
       do i=1,size(App_full,1)
          App_full(i,i)=-1.0_mp
       end do
    end if
    call Full2Sparse(App_full,problem%App)
  end subroutine init_acoustic_operator


  ! Initializes the main operators with ABC
   subroutine init_acoustic_operator_abc(problem)
    type(acoustic_problem),intent(inout)        :: problem
    real(mp),dimension(problem%nb_elem*problem%DoF,                             &
                   problem%nb_elem*problem%DoF) :: Ap_full,Av_full,B,App_full,  &
                                                   m_glob,minv_glob,s_glob,     &
                                                   Fp,Fv,Dp,Dv,dummy,           &
                                                   minv_glob_abc,mabs
    real(mp),dimension(problem%DoF,problem%DoF) :: m_loc,minv_loc,s_loc
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
    call init_FpFv(Fp,Fv,DoF,nb_elem,problem%boundaries,problem%flux,           &
                   problem%alpha,problem%receiver_loc)
    call init_DpDv(Dp,Dv,DoF,nb_elem,problem%velocity,problem%density)

    Ap_full=0.0_mp
    Av_full=0.0_mp
    App_full=0.0_mp
    call Full2Sparse(m_glob,problem%M)
    call Full2Sparse(minv_glob,problem%Minv)
    call Full2Sparse(s_glob,problem%S)

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
          Av_full=(1.0_mp/24)*(problem%dt/problem%dx)**2*matmul(Av_full,B)+Av_full
          Ap_full=(1.0_mp/24)*(problem%dt/problem%dx)**2*matmul(B,Ap_full)+Ap_full
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
          Av_full=(1.0_mp/24)*(problem%dt/problem%dx)**2*matmul(Av_full,B)+Av_full
          Ap_full=(1.0_mp/24)*(problem%dt/problem%dx)**2*matmul(B,Ap_full)+Ap_full
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
  end subroutine init_acoustic_operator_abc


  !**************** PROBLEM RESOLUTION **** *************************************
  ! Makes one forward time step
  subroutine one_forward_step(problem,t)
    type(acoustic_problem),intent(inout) :: problem
    real(mp)              ,intent(in)    :: t

    real(mp),dimension(problem%DoF*problem%nb_elem) :: RHSu0,RHSp0
    real(mp),dimension(problem%DoF*problem%nb_elem) :: RHSu1,RHSp1
    real(mp),dimension(problem%DoF*problem%nb_elem) :: RHSu2,RHSp2
    real(mp),dimension(problem%DoF*problem%nb_elem) :: U1,U2,U3,U4
    real(mp),dimension(problem%DoF*problem%nb_elem) :: P1,P2,P3,P4
    real(mp)                                        :: gd,gg
    real(mp)                                        :: f0
    integer                                         :: i
    real(mp),dimension(problem%DoF*problem%nb_elem) :: FP_0,FP_half,FP_1
    real(mp),dimension(problem%DoF*problem%nb_elem) :: FU_0,FU_half,FU_1

    if ((problem%time_scheme.eq.'LF').or.(problem%time_scheme.eq.'LF4')) then

       call eval_RHS(FP_0,FU_0,problem,t-problem%dt)
       call eval_RHS(FP_half,FU_half,problem,t-0.5_mp*problem%dt)

       call LF_forward(problem%P,problem%U,problem%Ap,problem%Av,problem%App,   &
            FP_0,FU_half)

    else if (problem%time_scheme.eq.'RK4') then

       
       call eval_RHS(FP_0,FU_0,problem,t-problem%dt)

       call eval_RHS(FP_half,FU_half,problem,t-0.5_mp*problem%dt)
       problem%RHSu_half=FU_half 
       problem%RHSp_half=FP_half      
       call eval_RHS(FP_1,FU_1,problem,t)
       problem%RHSu=FU_1
       problem%RHSp=FP_1

       call RK4_forward(problem%P,problem%U,problem%Ap,problem%Av,problem%App,  &
                        FP_0,FP_half,FP_1,FU_0,FU_half,FU_1)

    else if (problem%time_scheme.eq.'AB3') then

       RHSp0=0.0_mp
       RHSu0=0.0_mp
       RHSp1=0.0_mp
       RHSu1=0.0_mp
       RHSp2=0.0_mp
       RHSu2=0.0_mp

       call eval_RHS(RHSp0,RHSu0,problem,t-problem%dt)
       problem%RHSu=RHSu0
       problem%RHSp=RHSp0
       call eval_RHS(RHSp2,RHSu1,problem,t-2*problem%dt)
       call eval_RHS(RHSp2,RHSu2,problem,t-3*problem%dt)
       
       call AB3_forward(problem%P,problem%U,problem%Ap,problem%Av,problem%App,  &
                      RHSp0,RHSu0,                                              &
                      RHSp1,RHSu1,                                              &
                      RHSp2,RHSu2,                                              &
                      problem%Pk1,problem%Pk2,problem%Uk1,problem%Uk2)

    end if
  end subroutine one_forward_step

  ! Makes all the forward steps
  subroutine all_forward_step(problem)
    type(acoustic_problem),intent(inout) :: problem
    integer                              :: i
    real(mp)                             :: t

    do i=1,problem%n_time_step
       t=i*problem%dt
       call one_forward_step(problem,t)
    end do
  end subroutine all_forward_step


  ! Frees the acoustic problem variables
  subroutine free_acoustic_problem(problem)
    type(acoustic_problem),intent(inout) :: problem
    deallocate(problem%U)
    deallocate(problem%P)
    deallocate(problem%density)
    deallocate(problem%velocity)
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
    call free_sparse_matrix(problem%S)
    call free_sparse_matrix(problem%App)
  end subroutine free_acoustic_problem


  ! Evaluates the RHS for the forward problem at time=t
  subroutine eval_RHS(RHSp,RHSu,problem,t)
    real(mp),dimension(:)   ,intent(inout) :: RHSp
    real(mp),dimension(:)   ,intent(inout) :: RHSu
    type(acoustic_problem)  ,intent(in)    :: problem
    real(mp)                ,intent(in)    :: t

    real(mp)    :: gg,gd,f0
    integer :: last_node

    RHSp=0.0_mp
    RHSu=0.0_mp
    f0=2.5_mp

    last_node=problem%DoF*problem%nb_elem

    if (t.le.0.0_mp) then
       RHSp=0.0_mp
       RHSu=0.0_mp
    else

       if (problem%boundaries.eq.'dirichlet') then
          gg=min(0.5_mp*t,0.5_mp)
          gd=-gg
          RHSu(1)=gg*(-1.0_mp+0.0_mp*problem%alpha)
          RHSu(last_node)=gd*(1.0_mp-0.0_mp*problem%alpha)
       end if

       if (problem%signal.eq.'flat') then
          gg=(2*t-1/f0)*exp(-(2.0_mp*PI*(2*t-1/f0)*f0)**2.0_mp)*5_mp/0.341238111_mp
          !gg=sin(4*PI*t)
          !gg=1.0_mp
          RHSp((problem%source_loc-1)*problem%DoF+1)=gg*(1.0_mp-0.0_mp*problem%alpha)
       end if
    end if

    RHSp=sparse_matmul(problem%Minv_p,RHSp)
    RHSu=sparse_matmul(problem%Minv_v,RHSu)
  end subroutine eval_RHS



  ! Evaluates the Error for an homogeneous periodic test case 
  subroutine error_periodic(problem,t,errorU,errorP,iter)
    type(acoustic_problem)    ,intent(in)  :: problem
    real(mp)                  ,intent(in)  :: t
    real(mp)                  ,intent(out) :: errorU
    real(mp)                  ,intent(out) :: errorP
    integer,optional          ,intent(in)  :: iter

    real(mp),dimension(problem%nb_elem*problem%DoF) :: U_ex
    real(mp),dimension(problem%nb_elem*problem%DoF) :: U_work
    real(mp),dimension(problem%nb_elem*problem%DoF) :: P_ex
    real(mp),dimension(problem%nb_elem*problem%DoF) :: P_work
    real(mp)                                        :: dx,ddx,x
    integer                                         :: i,j,k
    integer                                         :: DoF,nb_elem

    if (problem%boundaries.ne.'periodic') then
       print*,'The error cannot be calculated'
       print*,'The boundary conditions are not periodic'
       RETURN
    end if

    DoF=problem%DoF
    nb_elem=problem%nb_elem

    dx=problem%total_length/nb_elem
    ddx=dx/(DoF-1)

    do i=1,nb_elem
       do j=1,DoF
          x=(i-1)*problem%dx+(j-1)*ddx-t

          do while (x.le.0.0_mp)
             x=x+problem%total_length
          end do

          U_ex((i-1)*DoF+j)=initial_perturbation(x,0.0_mp,problem%signal)
          P_ex((i-1)*DoF+j)=initial_perturbation(x,-problem%dt,problem%signal)
       end do
    end do

    errorU=0.0_mp
    errorP=0.0_mp

    if(problem%bernstein) then
       do k=1,nb_elem
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

    errorU=sqrt(errorU)/sqrt(dot_product(U_ex,U_ex))
    errorP=sqrt(errorP)/sqrt(dot_product(P_ex,P_ex))

    print*,'-----------------------------'
    print*,'-----------------------------'
    print*,'t=',t
    print*,'Lerreur en U est de : ', errorU
    print*,'Lerreur en P est de : ', errorP
    print*,'-----------------------------'    
    print*,'-----------------------------'
  end subroutine error_periodic  
end module m_acoustic
