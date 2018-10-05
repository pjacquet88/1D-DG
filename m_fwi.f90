module m_fwi
  use m_file_function
  use m_matrix
  use m_time_scheme
  use m_acoustic
  use m_adjoint
  implicit none


  type t_fwi
     type(acoustic_problem)           :: forward
     type(adjoint_problem)            :: backward
     integer                          :: nb_iter
     integer                          :: model_size
     real,dimension(:)  ,allocatable  :: gradJ_velocity
     real,dimension(:)  ,allocatable  :: gradJ_density
     real,dimension(:)  ,allocatable  :: velocity_model
     real,dimension(:)  ,allocatable  :: density_model
     integer                          :: n_time_step
     real,dimension(:,:),allocatable  :: data_P
     real,dimension(:,:),allocatable  :: data_U
     real,dimension(:,:),allocatable  :: P_received
     real,dimension(:,:),allocatable  :: U_received
     real,dimension(:,:),allocatable  :: P
     real,dimension(:,:),allocatable  :: U
     real,dimension(:,:),allocatable  :: QP
     real,dimension(:,:),allocatable  :: QU
     
     integer                          :: nb_elem
     integer                          :: DoF
     character(len=20)                :: time_scheme
     real                             :: total_length
     real                             :: final_time
     real                             :: alpha
     logical                          :: bernstein
     character(len=20)                :: signal
     character(len=20)                :: boundaries
     integer                          :: k_max
     real                             :: epsilon
     integer                          :: source_loc
     integer                          :: receiver_loc
  end type t_fwi
    

contains

  subroutine init_fwi(fwi,nb_iter,velocity_model,density_model,data_P,data_U,   &
                      nb_elem,DoF,time_scheme,total_length,final_time,alpha,    &
                      bernstein,k_max,epsilon,source_loc,receiver_loc)
    type(t_fwi)        ,intent(inout) :: fwi
    integer            ,intent(in)    :: nb_iter
    real,dimension(:)  ,intent(in)    :: velocity_model
    real,dimension(:)  ,intent(in)    :: density_model
    real,dimension(:,:),intent(in)    :: data_P
    real,dimension(:,:),intent(in)    :: data_U
    integer            ,intent(in)    :: nb_elem
    integer            ,intent(in)    :: DoF
    character(len=*)   ,intent(in)    :: time_scheme
    real               ,intent(in)    :: total_length
    real               ,intent(in)    :: final_time
    real               ,intent(in)    :: alpha
    logical            ,intent(in)    :: bernstein
    integer            ,intent(in)    :: k_max
    real               ,intent(in)    :: epsilon
    integer            ,intent(in)    :: source_loc
    integer            ,intent(in)    :: receiver_loc
    
    fwi%nb_iter=nb_iter
    fwi%model_size=size(density_model)

    allocate(fwi%velocity_model(fwi%model_size))
    allocate(fwi%density_model(fwi%model_size))

    fwi%velocity_model=velocity_model
    fwi%density_model=density_model

    allocate(fwi%data_P(0:size(data_P,1)-1,size(data_U,2)))
    allocate(fwi%data_U(0:size(data_U,1)-1,size(data_U,2)))

    fwi%data_P=data_P
    fwi%data_U=data_U

    fwi%nb_elem=nb_elem
    fwi%DoF=DoF
    fwi%time_scheme=time_scheme
    fwi%total_length=total_length
    fwi%final_time=final_time
    fwi%alpha=alpha
    fwi%bernstein=bernstein
    fwi%k_max=k_max
    fwi%epsilon=epsilon
    fwi%source_loc=source_loc
    fwi%receiver_loc=receiver_loc
    
    allocate(fwi%gradJ_velocity(fwi%model_size))
    allocate(fwi%gradJ_density(fwi%model_size))

    fwi%signal='flat'
    fwi%boundaries='ABC'
  end subroutine init_fwi
  

  subroutine free_fwi(fwi)
    type(t_fwi),intent(inout) :: fwi

    deallocate(fwi%velocity_model)
    deallocate(fwi%density_model)
    deallocate(fwi%data_P)
    deallocate(fwi%data_U)
    deallocate(fwi%gradJ_velocity)
    deallocate(fwi%gradJ_density)
  end subroutine free_fwi


  subroutine one_fwi_step(fwi)
    type(t_fwi),intent(inout) :: fwi
    integer                   :: current_time_step,i
    real                      :: t

    !------------------------------ Forward -------------------------------------
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%% FORWARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call init_acoustic_problem(fwi%forward,fwi%nb_elem,fwi%DoF,fwi%time_scheme, &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%k_max,           &
         fwi%epsilon,fwi%source_loc,fwi%receiver_loc)



    print*,'Initialization done'
    
    
    fwi%n_time_step=fwi%forward%n_time_step


    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%% backward %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call init_adjoint_problem(fwi%backward,fwi%nb_elem,fwi%DoF,fwi%time_scheme, &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%k_max,           &
         fwi%epsilon,fwi%source_loc,fwi%receiver_loc)

        print*,'Il y a :',fwi%n_time_step,' time steps'
    
    allocate(fwi%P(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%U(fwi%DoF*fwi%nb_elem,0:fwi%n_time_step))
    allocate(fwi%P_received(0:fwi%n_time_step,2))
    allocate(fwi%U_received(0:fwi%n_time_step,2))
    fwi%P(:,0)=0
    fwi%U(:,0)=0

    current_time_step=0
    t=0
    do current_time_step=1,fwi%n_time_step
       t=current_time_step*fwi%forward%dt
       call one_forward_step(fwi%forward,t)
       fwi%P(:,current_time_step)=fwi%forward%P
       fwi%U(:,current_time_step)=fwi%forward%U
       fwi%P_received(current_time_step,1)=t
       fwi%P_received(current_time_step,2)=fwi%forward%P(fwi%receiver_loc)
       fwi%U_received(current_time_step,1)=t
       fwi%U_received(current_time_step,2)=fwi%forward%P(fwi%receiver_loc)
    end do

    print*,'Calculus done'

    do i=0,fwi%n_time_step
       if (modulo(i,100).eq.0) then
          call print_vect(fwi%P(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'P')
       else if (i.eq.0) then
          call print_vect(fwi%P(:,i),fwi%nb_elem,fwi%DoF,fwi%forward%dx,fwi%bernstein,i,'P')
       end if
    end do



    call free_acoustic_problem(fwi%forward)

  end subroutine one_fwi_step


end module m_fwi
