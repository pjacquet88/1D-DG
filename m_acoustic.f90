module m_acoustic
  use m_polynom
  use m_matrix
  use m_powermethod
  use m_file_function
  implicit none
  
  type acoustic_problem
     integer                         :: DoF          ! nb of degree of freedom 
     integer                         :: nb_elem      ! nb of elements
     integer                         :: time_order   ! time order 2 or 4 (LF)
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
     character(len=20)               :: boundaries   ! string for BC
     character(len=20)               :: signal       ! string to initialize U,P
     integer                         :: k_max        ! iter max for power method
     real                            :: epsilon      ! precition for power method algo.
     !----------- MATRICES ---------------------------
     type(sparse_matrix)             :: Ap,App       
     type(sparse_matrix)             :: Av
     type(sparse_matrix)             :: Minv_p,Minv_v
     !----------- SOURCE & RECEIVER ---------------------------
     integer                         :: source_loc   ! beginning of an element
     integer                         :: receiver_loc ! beginning of an element
     real,dimension(:,:),allocatable :: receiver     ! vector describing receivers values as a function of time
     real,dimension(:,:),allocatable :: data         ! vector describing real datas values as a function of time
     real,dimension(:,:),allocatable :: RHS_backward ! what we inject in the backward
     integer                         :: n_time_step  ! nb of time step
     integer                         :: n_display    ! nb of time step between each saves
     integer                         :: n_frame      ! nb of time where the sol. is saved
     logical                         :: forward      ! T=forward of F=backward
  end type acoustic_problem
  
  real,dimension(15)                 :: xx,weight
  real,parameter :: PI=acos(-1.0)

  public  :: init_problem,init_operator,all_time_step,                          &
             free_problem,                                                      &
             print_sol,error_periodique

  private :: xx,weight,solution,PI,                                             &
             init_quadrature,Int_bxb,Int_lxl,init_m_loc,init_s_loc,             &
             init_UP,init_m_glob,init_minv_glob,init_minv_glob_abc,init_mabs,   &
             init_s_glob,init_FpFv,init_DpDv,init_dt,eval_RHS,one_time_step
  
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
  
  function solution(x,dt,signal)
    real,intent(in) :: x
    real,intent(in) :: dt
    character(len=*),intent(in) :: signal
    real                        :: solution

    real                        :: xt

    xt=x+dt
    
    if (signal.eq.'sinus') then
       solution=sin(4*PI*xt)
       ! solution=sin(xt)
    else if (signal.eq.'creneau') then
       if ((x.lt.0.55).and.(x.gt.0.45)) then
          solution=0.5
       else
          solution=0.0
       end if
    else if (signal.eq.'x') then
       solution=xt**5
        else if (signal.eq.'flat') then
       solution=0.0
    else
       solution=(xt-0.5)*exp(-(2.0*PI*(xt-0.5)*2.0)**2.0)*20
    end if    
  end function solution


  subroutine init_problem(problem,nb_elem,DoF,time_order,velocity,density,      &
                          total_length,final_time,alpha,bernstein,signal,       &
                          boundaries,k_max,epsilon,source_loc,receiver_loc,     &
                          n_frame,forward)
    
    type(acoustic_problem),intent(out) :: problem
    integer               ,intent(in)  :: nb_elem
    integer               ,intent(in)  :: DoF
    real,dimension(:)     ,intent(in)  :: velocity
    real,dimension(:)     ,intent(in)  :: density
    integer               ,intent(in)  :: time_order
    real                  ,intent(in)  :: total_length
    real                  ,intent(in)  :: final_time
    real                  ,intent(in)  :: alpha
    logical               ,intent(in)  :: bernstein
    character(len=*)      ,intent(in)  :: signal
    character(len=*)      ,intent(in)  :: boundaries
    integer               ,intent(in)  :: k_max
    real                  ,intent(in)  :: epsilon
    integer               ,intent(in)  :: source_loc
    integer               ,intent(in)  :: receiver_loc
    integer               ,intent(in)  :: n_frame
    logical               ,intent(in)  :: forward
    
    integer :: i,j
    integer :: dv,dp,size_velocity,size_density
    real    :: gd,dd

    problem%DoF=DoF
    problem%time_order=time_order
    problem%nb_elem=nb_elem
    problem%total_length=total_length
    problem%dx=total_length/(nb_elem)
    problem%final_time=final_time
    problem%alpha=alpha
    problem%bernstein=bernstein
    problem%boundaries=boundaries
    problem%k_max=k_max
    problem%epsilon=epsilon
    problem%signal=signal
    problem%receiver_loc=receiver_loc
    problem%source_loc=source_loc
    problem%n_frame=n_frame
    problem%forward=forward

    
    size_velocity=size(velocity)
    size_density=size(density)
    
    allocate(problem%U(DoF*nb_elem),problem%P(DoF*nb_elem))
    problem%U=0.0
    problem%P=0.0
    if (size(velocity).gt.nb_elem) then
       print*,'Velocity input size does not correspond'
       STOP
    end if
    if (size(density).gt.nb_elem) then
       print*,'Density input size does not correspond'
       STOP
    end if
    allocate(problem%velocity(DoF*nb_elem),problem%density(DoF*nb_elem))

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
       call init_operator_abc(problem)
    else
       call init_operator(problem)
    end if
    
    problem%n_time_step=int(problem%final_time/problem%dt)
    if (problem%n_time_step.le.problem%n_frame) then
       problem%n_display=1 ! to avoid errors
    else
       problem%n_display=int(problem%n_time_step/problem%n_frame)+1
    end if

    problem%n_time_step=problem%n_display*problem%n_frame
    problem%final_time=problem%dt*problem%n_time_step
    ! here I change the final time to fit frames of the forward and the backward
    
    allocate(problem%receiver(0:problem%n_time_step/problem%n_display,2))
    if (.not.problem%forward) then
       allocate(problem%data(0:problem%n_time_step/problem%n_display,2))
       allocate(problem%RHS_backward(0:problem%n_time_step/problem%n_display,2))
       open(unit=10,file='data.dat',status="old", action="read")
       do i=0,size(problem%data,1)-1
          read(10,*) problem%data(i,1),problem%data(i,2)
       end do
       close(10)
       open(unit=10,file='receiver.dat',status="old", action="read")
       do i=0,size(problem%RHS_backward,1)-1
          read(10,*) problem%RHS_backward(i,1),problem%RHS_backward(i,2)
       end do
       close(10)
       do i=0,size(problem%RHS_backward,1)-1
          problem%RHS_backward(i,2)=problem%data(i,2)*                          &
               min(1.0,10*exp(-50*(problem%data(i,1)/problem%final_time-0.8)*   &
               problem%data(i,1)/problem%final_time))!-problem%RHS_backward(i,2)
       end do
       open(unit=35,file='data_backward.dat')
       do i=0,size(problem%RHS_backward,1)-1
          write(35,*) problem%RHS_backward(i,1),problem%RHS_backward(i,2)
       end do
    end if

    print*,'Ap Av done'
    call init_UP(problem)
    print*,'nb_elem              ::',nb_elem
    print*,'DoF                  ::',DoF
    print*,'total_length         ::',total_length
    print*,'final_time           ::',problem%final_time
    print*,'dx                   ::',problem%dx
    print*,'dt                   ::',problem%dt
    print*,'n_time_step          ::',problem%n_time_step
    print*,'n_time_step          ::',problem%n_display
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
    print*,'------------------------------------------------------------'    
  end subroutine init_problem

  
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
  

 !--------------------- Init Matrices ------------------------------------

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

  
  subroutine init_minv_loc(minv_loc,m_loc)
    real,dimension(:,:),intent(out) :: minv_loc
    real,dimension(:,:),intent(in)  :: m_loc
    minv_loc=LU_inv(m_loc)
  end subroutine init_minv_loc

  
  subroutine init_m_glob(m_glob,m_loc,nb_elem)
    real,dimension(:,:),intent(out) :: m_glob
    real,dimension(:,:),intent(in)  :: m_loc
    integer            ,intent(in)  :: nb_elem

    integer :: i,DoF

    DoF=size(m_loc,1)
        
    do i=1,nb_elem
       m_glob(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=m_loc
    end do
  end subroutine init_m_glob
  
  
  subroutine init_minv_glob(minv_glob,minv_loc,nb_elem)
    real,dimension(:,:),intent(out) :: minv_glob       
    real,dimension(:,:),intent(in)  :: minv_loc
    integer            ,intent(in)  :: nb_elem
    integer                         :: i,DoF

    DoF=size(minv_loc,1)
    
    do i=1,nb_elem
       minv_glob(DoF*(i-1)+1:DoF*i,DoF*(i-1)+1:DoF*i)=minv_loc
    end do
    
  end subroutine init_minv_glob

  subroutine init_minv_glob_abc(minv_glob_abc,m_loc,nb_elem,DoF,velocity_beg,   &
                                velocity_end,dx,dt)
    real,dimension(:,:),intent(out) :: minv_glob_abc
    real,dimension(:,:),intent(in)  :: m_loc
    integer            ,intent(in)  :: nb_elem
    integer            ,intent(in)  :: DoF
    real               ,intent(in)  :: velocity_beg
    real               ,intent(in)  :: velocity_end
    real               ,intent(in)  :: dx
    real               ,intent(in)  :: dt

    integer                         :: i
    real,dimension(DoF,DoF)         :: m_loc_abc1,m_loc_abc2
    real,dimension(DoF,DoF)         :: minv_loc,minv_loc_abc1,minv_loc_abc2

    m_loc_abc1=m_loc
    m_loc_abc2=m_loc
    m_loc_abc2(DoF,DoF)=m_loc_abc2(DoF,DoF)+0.5*velocity_end*dt/dx
    m_loc_abc1(1,1)=m_loc_abc1(1,1)+0.5*velocity_beg*dt/dx

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

  subroutine init_mabs(mabs,nb_elem,DoF,velocity_beg,velocity_end,dx,dt)
    real,dimension(:,:),intent(out) :: mabs
    integer            ,intent(in)  :: nb_elem
    integer            ,intent(in)  :: DoF
    real               ,intent(in)  :: velocity_beg
    real               ,intent(in)  :: velocity_end
    real               ,intent(in)  :: dx
    real               ,intent(in)  :: dt

    integer :: i
    
    mabs=0.0
    mabs(DoF*nb_elem,DoF*nb_elem)=0.5*velocity_end*dt/dx
    mabs(1,1)=0.5*velocity_beg*dt/dx
  end subroutine init_mabs

    
  subroutine init_s_loc(s_loc,DoF,bernstein)
    real,dimension(:,:),intent(out) :: s_loc
    integer            ,intent(in)  :: DoF
    logical            ,intent(in)  :: bernstein
    integer                         :: i,j
    real,dimension(DoF)             :: bpol,dbpol
    real,dimension(DoF)             :: dlpol

    
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

  
  subroutine init_s_glob(s_glob,s_loc,nb_elem)
    real,dimension(:,:),intent(out) :: s_glob
    real,dimension(:,:),intent(in)  :: s_loc
    integer            ,intent(in)  :: nb_elem
    integer :: i,DoF

    DoF=size(s_loc,1)


    do i=1,nb_elem
       s_glob(DoF*(i-1)+1:Dof*i,DoF*(i-1)+1:Dof*i)=s_loc
    end do
  end subroutine init_s_glob

  subroutine init_FpFv(Fp,Fv,DoF,nb_elem,boundaries,alpha,forward,receiver_loc)
    real,dimension(:,:),intent(out) :: Fp
    real,dimension(:,:),intent(out) :: Fv
    integer            ,intent(in)  :: DoF
    integer            ,intent(in)  :: nb_elem
    character(len=*)   ,intent(in)  :: boundaries
    real               ,intent(in)  :: alpha
    logical            ,intent(in)  :: forward
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

    if (boundaries.eq.'periodique') then
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

       Fp(1,1)=1.0+2.0*alpha_dummy
       Fp(1,DoF*nb_elem)=0.0

       Fp(Dof*nb_elem,Dof*nb_elem)=-1.0-2.0*alpha_dummy
       Fp(Dof*nb_elem,1)=0.0

       alpha_dummy=-alpha

       Fv(1,1)=0.0
       Fv(1,DoF*nb_elem)=0.0

       Fv(Dof*nb_elem,Dof*nb_elem)=0.0
       Fv(Dof*nb_elem,1)=0.0

    elseif (boundaries.eq.'ABC') then

       alpha_dummy=alpha

       Fv(1,1)=0.0!(1.0+2.0*alpha_dummy)
       Fv(1,DoF*nb_elem)=0.0

       Fv(Dof*nb_elem,Dof*nb_elem)=0.0!-1.0-2.0*alpha_dummy
       Fv(Dof*nb_elem,1)=0.0

       ! if (.not.forward) then
       !    Fv(DoF*receiver_loc,DoF*receiver_loc)=-1.0-2.0*alpha_dummy
       !    Fv(DoF*receiver_loc+1,DoF*receiver_loc+1)=1.0+2.0*alpha_dummy
       ! end if
       
       alpha_dummy=-alpha

!       Fp(1,1)=1.0
       Fp(1,1)=1.0
       Fp(1,DoF*nb_elem)=0.0

       Fp(Dof*nb_elem,Dof*nb_elem)=-1.0!*2
       Fp(Dof*nb_elem,1)=0.0

    elseif (boundaries.eq.'neumann') then
       alpha_dummy=alpha

       Fp(1,1)=0.0
       Fp(1,DoF*nb_elem)=0.0

       Fp(Dof*nb_elem,Dof*nb_elem)=0.0
       Fp(Dof*nb_elem,1)=0.0

       alpha_dummy=-alpha

       Fv(1,1)=1.0-2.0*alpha_dummy
       Fv(1,DoF*nb_elem)=0.0
       
       Fv(Dof*nb_elem,Dof*nb_elem)=-1.0+2.0*alpha_dummy
       Fv(Dof*nb_elem,1)=0.0
    end if
  end subroutine init_FpFv


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


  subroutine init_dt(Ap,Av,dt,k_max,epsilon,nb_elem,time_order)
    real,dimension(:,:),intent(in)  :: Ap
    real,dimension(:,:),intent(in)  :: Av
    real               ,intent(out) :: dt
    integer            ,intent(in)  :: k_max
    real               ,intent(in)  :: epsilon
    integer            ,intent(in)  :: nb_elem
    integer            ,intent(in)  :: time_order
    type(sparse_matrix)             :: sparse_A
    real,dimension(size(Ap,1),size(Ap,2)) :: A
    real                            :: max_value
    integer                         :: i
    real,parameter                  :: alpha=1.00
    real                            :: cfl

    A=matmul(Ap,Av)
    ! do i=1,size(A,1)
    !    write(789,*) A(i,:)
    ! end do
    
    call Full2Sparse(A,sparse_A)
    call power_method_sparse(sparse_A,max_value,k_max,epsilon)
    dt=alpha/(sqrt(abs(max_value)))
    call free_sparse_matrix(sparse_A)
    A=matmul(Av,Ap)
    call Full2Sparse(matmul(Av,Ap),sparse_A)
    call power_method_sparse(sparse_A,max_value,k_max,epsilon)
    dt=min(dt,alpha/(sqrt(abs(max_value))))
    call free_sparse_matrix(sparse_A)
    dt=dt/2.0

    !****** dt constant *************************
    !****** 100 ***************************
    ! dt= 4.16666589E-04*alpha  !ordre  1 
    ! dt= 2.08152022E-04*alpha  !ordre  2
    ! dt= 1.24648141E-04*alpha  !ordre  3
    ! dt= 8.29414494E-05*alpha  !ordre  4
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
     ! cfl=0.95*dt*100 ! cfl=dt/h
     ! if (time_order.eq.4) then
     !    cfl=cfl*2.7
     ! end if
     ! dt=cfl/(nb_elem)
  end subroutine init_dt
  

  subroutine init_operator(problem)
    type(acoustic_problem),intent(inout)        :: problem
    real,dimension(problem%nb_elem*problem%DoF,                                 &
                   problem%nb_elem*problem%DoF) :: Ap_full,Av_full,B,App_full,  &
                                                   m_glob,minv_glob,s_glob,     &
                                                   Fp,Fv,Dp,Dv
    real,dimension(problem%DoF,problem%DoF)     :: m_loc,minv_loc,s_loc
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
    call init_FpFv(Fp,Fv,DoF,nb_elem,problem%boundaries,problem%alpha,          &
                   problem%forward,problem%receiver_loc)
    call init_DpDv(Dp,Dv,DoF,nb_elem,problem%velocity,problem%density)
    
    Ap_full=0.0
    Av_full=0.0

    if (problem%bernstein) then

       Ap_full=(s_glob+matmul(minv_glob,Fp))
       Av_full=(s_glob+matmul(minv_glob,Fv))
       Ap_full=matmul(Dp,Ap_full)
       Av_full=matmul(Dv,Av_full)

       call init_dt((1/problem%dx)*Ap_full,(1/problem%dx)*Av_full,problem%dt,   &
            problem%k_max,problem%epsilon,problem%nb_elem,             &
            problem%time_order)

       if (problem%time_order.eq.4) then
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
                     problem%k_max,problem%epsilon,problem%nb_elem,             &
                     problem%time_order)

       if (problem%time_order.eq.4) then
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
    
    print*,'Non-Zero Values of Ap,Av  ::',problem%Ap%NNN
    print*,'Size of Ap,AV matrices    ::',problem%Ap%nb_ligne,'x', &
         problem%Ap%nb_ligne,'=',problem%Ap%nb_ligne**2
    print*,'Ratio                     ::',real(problem%Ap%NNN)/                 &
                                                           problem%Ap%nb_ligne**2
  end subroutine init_operator


  
   subroutine init_operator_abc(problem)
    type(acoustic_problem),intent(inout)        :: problem
    real,dimension(problem%nb_elem*problem%DoF,                                 &
                   problem%nb_elem*problem%DoF) :: Ap_full,Av_full,B,App_full,  &
                                                   m_glob,minv_glob,s_glob,     &
                                                   Fp,Fv,Dp,Dv,                 &
                                                   minv_glob_abc,mabs
    real,dimension(problem%DoF,problem%DoF)     :: m_loc,minv_loc,s_loc
    integer                                     :: i,j
    integer                                     :: DoF,nb_elem
    
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
                   problem%forward,problem%receiver_loc)
    call init_DpDv(Dp,Dv,DoF,nb_elem,problem%velocity,problem%density)

    Ap_full=0.0
    Av_full=0.0
    App_full=0.0
    
    if (problem%bernstein) then

       Av_full=(s_glob+matmul(minv_glob,Fv))
       Av_full=matmul(Dv,Av_full)

       call init_dt((1/problem%dx)*Av_full,(1/problem%dx)*Av_full,problem%dt,   &
            problem%k_max,problem%epsilon,problem%nb_elem,problem%time_order)

       call init_minv_glob_abc(minv_glob_abc,m_loc,nb_elem,DoF,                 &
                               problem%velocity(1),problem%velocity(nb_elem),   &
                               problem%dx,problem%dt)
       call init_mabs(mabs,nb_elem,DoF,                                         &
                      problem%velocity(1),problem%velocity(nb_elem),            &
                      problem%dx,problem%dt)
       
       Ap_full=matmul(m_glob,s_glob)+Fp
       Ap_full=matmul(minv_glob_abc,Ap_full)
       Ap_full=matmul(Dp,Ap_full)
      
       App_full=mabs-m_glob
       App_full=matmul(minv_glob_abc,App_full)

       if (problem%time_order.eq.4) then
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
                     problem%k_max,problem%epsilon,problem%nb_elem,             &
                     problem%time_order)


       call init_minv_glob_abc(minv_glob_abc,m_loc,nb_elem,DoF,                 &
                               problem%velocity(nb_elem),                       &
                               problem%density(nb_elem),problem%dx,problem%dt)
       call init_mabs(mabs,nb_elem,DoF,                         &
                              problem%velocity(nb_elem),                        &
                              problem%density(nb_elem),problem%dx,problem%dt)
       
       Ap_full=s_glob+Fp
       
       Ap_full=matmul(minv_glob_abc,Ap_full)
       Ap_full=matmul(Dp,Ap_full)
       
       
       App_full=mabs-m_glob
       App_full=matmul(minv_glob_abc,App_full)

       if (problem%time_order.eq.4) then
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
    
    print*,'Non-Zero Values of Ap,Av  ::',problem%Ap%NNN
    print*,'Size of Ap,AV matrices    ::',problem%Ap%nb_ligne,'x', &
         problem%Ap%nb_ligne,'=',problem%Ap%nb_ligne**2
    print*,'Ratio                     ::',real(problem%Ap%NNN)/problem%Ap%nb_ligne**2
  end subroutine init_operator_abc
  
  
 !**************** RESOLUTION DU PROBLEM *************************************
  subroutine one_time_step(problem,t)
    type(acoustic_problem),intent(inout) :: problem
    real,                  intent(in)    :: t

    real,dimension(problem%DoF*problem%nb_elem) :: RHSv,RHSp
    real                                        :: gd,gg
    real                                        :: f0
    integer                                     :: last_node,i

    last_node=problem%DoF*problem%nb_elem
    RHSp=0.0
    RHSv=0.0
    f0=2.5

    if (problem%boundaries.eq.'dirichlet') then
       gg=min(0.5*t,0.5)
       gd=0.0!-gg
       RHSp(1)=gg*(-1.0-2.0*problem%alpha)
       RHSp(last_node)=gd*(1.0+2.0*problem%alpha)
       
       RHSv(1)=0.0
       RHSv(size(RHSv))=0.0

       RHSv=sparse_matmul(problem%Minv_v,RHSv)
       RHSp=sparse_matmul(problem%Minv_p,RHSp)

    end if
    if (problem%boundaries.eq.'ABC') then
       if (problem%forward) then
          gg=(2*t-1/f0)*exp(-(2.0*PI*(2*t-1/f0)*f0)**2.0)*5/0.341238111
          RHSv((problem%source_loc-1)*problem%DoF+1)=gg*(-1.0-0.0*problem%alpha)
         ! RHSv(problem%receiver_loc*problem%DoF)=gg*(-1.0+2.0*problem%alpha)
       else
          gg=eval_RHS(problem%RHS_backward,t,problem%final_time)
          ! RHSv(problem%receiver_loc*problem%DoF+1)=gg*(-1.0-2.0*problem%alpha)
          RHSv(problem%receiver_loc*problem%DoF)=gg*(-1.0-2.0*problem%alpha)
       end if
          
          RHSv=sparse_matmul(problem%Minv_v,RHSv)
          RHSp=sparse_matmul(problem%Minv_p,RHSp)

    end if
    
    problem%U=problem%U-sparse_matmul(problem%Av,problem%P)-RHSv
    if (problem%boundaries.eq.'ABC') then
       problem%P=-sparse_matmul(problem%Ap,problem%U)-                          &
                  sparse_matmul(problem%App,problem%P)-RHSp
    else
       problem%P=problem%P-sparse_matmul(problem%Ap,problem%U)-RHSp
    end if
  end subroutine one_time_step


  subroutine all_time_step(problem,sortie)
    type(acoustic_problem),intent(inout) :: problem
    logical               ,intent(in)    :: sortie
    integer                              :: i
    real                                 :: t
    
    if (sortie) then
       call print_sol(problem,0)
    end if

    problem%receiver(0,1)=0.0
    problem%receiver(0,2)=0.0
    
    do i=1,problem%n_time_step
       t=i*problem%dt
       call one_time_step(problem,t)

       if (modulo(i,problem%n_display).eq.0) then

          if (problem%receiver_loc.ne.0) then
              problem%receiver(i/problem%n_display,1)=t
              problem%receiver(i/problem%n_display,2)=problem%P(problem%receiver_loc*problem%DoF+1)

             if (sortie) then
                if (modulo(i,problem%n_display).eq.0) then
                   call print_sol(problem,i/problem%n_display)
                end if
             end if
          end if
       end if
    end do
    
    if (problem%forward) then
       open(23, file="receiver.dat", form="formatted")
       do i=0,size(problem%receiver,1)-1
          write(23,*) problem%receiver(i,1),problem%receiver(i,2)
       end do
       close(23)
    end if
    
    if (problem%forward) then
       do i=0,problem%n_time_step
          t=i*problem%dt
          write(25,*) t,eval_RHS(problem%receiver,t,problem%final_time)
       end do
    end if
  end subroutine all_time_step

  subroutine free_problem(problem)
    type(acoustic_problem),intent(inout) :: problem
    deallocate(problem%U)
    deallocate(problem%P)
    deallocate(problem%density)
    deallocate(problem%velocity)
    call free_sparse_matrix(problem%Ap)
    call free_sparse_matrix(problem%Av)
    call free_sparse_matrix(problem%Minv_p)
    call free_sparse_matrix(problem%Minv_v)

    if (problem%boundaries.eq.'ABC') then
       call free_sparse_matrix(problem%App)
    end if
  end subroutine free_problem

  subroutine print_sol(problem,N)
    type(acoustic_problem),intent(in) :: problem
    integer               ,intent(in) :: N

    if (problem%forward) then
      ! call print_vect(problem%U,problem%nb_elem,problem%DoF,problem%dx,           &
     !       problem%bernstein,N,'U')

       call print_vect(problem%P,problem%nb_elem,problem%DoF,problem%dx,           &
            problem%bernstein,N,'P')
    else
      ! call print_vect(problem%U,problem%nb_elem,problem%DoF,problem%dx,           &
      !      problem%bernstein,N,'Ub')

       call print_vect(problem%P,problem%nb_elem,problem%DoF,problem%dx,           &
            problem%bernstein,problem%n_frame-N,'B')
    end if
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

    if (problem%boundaries.ne.'periodique') then
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

          do while (x.le.0.0)
             x=x+problem%total_length
          end do
          
          U_ex((i-1)*DoF+j)=solution(x,0.0,problem%signal)
          P_ex((i-1)*DoF+j)=solution(x,-problem%dt,problem%signal)
       end do
    end do

    errorU=0
    errorP=0

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

  end subroutine error_periodique

  function eval_RHS(RHS,t,final_time)
    real,dimension(:,:),intent(in) :: RHS
    real               ,intent(in) :: t
    real               ,intent(in) :: final_time
    real                           :: eval_RHS
    integer                        :: i
    real                           :: t2

    t2=final_time-t

    i=1
    do while ((RHS(i,1).le.t2).and.(i.ne.size(RHS,1)))
       i=i+1
    end do
    i=i-1

    eval_RHS=(t2-RHS(i,1))*RHS(i+1,2)+(RHS(i+1,1)-t2)*RHS(i,2)
    eval_RHS=eval_RHS/(RHS(i+1,1)-RHS(i,1))

  end function eval_RHS
  
end module m_acoustic
