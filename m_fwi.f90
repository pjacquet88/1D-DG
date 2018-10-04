module m_fwi
  use m_file_function
  use m_matrix
  use m_time_scheme
  use m_acoustic
  implicit none


  type t_fwi
     type(acoustic_problem)           :: forward
     type(acoustic_problem)           :: backward
     integer                          :: nb_iter
     integer                          :: model_size
     real,dimension(:)  ,allocatable  :: gradJ_velocity
     real,dimension(:)  ,allocatable  :: gradJ_density
     real,dimension(:)  ,allocatable  :: velocity_model
     real,dimension(:)  ,allocatable  :: density_model
     integer                          :: nb_frame
     real,dimension(:)  ,allocatable  :: data
     real,dimension(:)  ,allocatable  :: P_received
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

  subroutine init_fwi(fwi,nb_iter,velocity_model,density_model,data,nb_elem,DoF,&
                      time_scheme,total_length,final_time,alpha,bernstein,k_max,&
                      epsilon,source_loc,receiver_loc)
    type(t_fwi)      ,intent(inout) :: fwi
    integer          ,intent(in)  :: nb_iter
    real,dimension(:),intent(in)  :: velocity_model
    real,dimension(:),intent(in)  :: density_model
    real,dimension(:),intent(in)  :: data
    integer          ,intent(in)  :: nb_elem
    integer          ,intent(in)  :: DoF
    character(len=*) ,intent(in)  :: time_scheme
    real             ,intent(in)  :: total_length
    real             ,intent(in)  :: final_time
    real             ,intent(in)  :: alpha
    logical          ,intent(in)  :: bernstein
    integer          ,intent(in)  :: k_max
    real             ,intent(in)  :: epsilon
    integer          ,intent(in)  :: source_loc
    integer          ,intent(in)  :: receiver_loc
    
    
    fwi%nb_iter=nb_iter
    fwi%model_size=size(density_model)

    allocate(fwi%velocity_model(fwi%model_size))
    allocate(fwi%density_model(fwi%model_size))

    fwi%velocity_model=velocity_model
    fwi%density_model=density_model

    fwi%nb_frame=size(data)
    allocate(fwi%data(fwi%nb_frame))

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
    allocate(fwi%P_received(fwi%nb_frame))
    allocate(fwi%P(fwi%model_size,fwi%nb_frame))
    allocate(fwi%U(fwi%model_size,fwi%nb_frame))
    allocate(fwi%QP(fwi%model_size,fwi%nb_frame))
    allocate(fwi%QU(fwi%model_size,fwi%nb_frame))

    fwi%signal='flat'
    fwi%boundaries='ABC'
  end subroutine init_fwi
  

  subroutine free_fwi(fwi)
    type(t_fwi),intent(inout) :: fwi

    deallocate(fwi%velocity_model)
    deallocate(fwi%density_model)
    deallocate(fwi%data)
    deallocate(fwi%gradJ_velocity)
    deallocate(fwi%gradJ_density)
    deallocate(fwi%P_received)
    deallocate(fwi%P)
    deallocate(fwi%U)
    deallocate(fwi%QP)
    deallocate(fwi%QU)
  end subroutine free_fwi


  subroutine one_fwi_step(fwi)
    type(t_fwi),intent(inout) :: fwi

    !------------------------------ Forward -------------------------------------
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%% FORWARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    call init_problem(fwi%forward,fwi%nb_elem,fwi%DoF,fwi%time_scheme,          &
         fwi%velocity_model,fwi%density_model,fwi%total_length,fwi%final_time,  &
         fwi%alpha,fwi%bernstein,fwi%signal,fwi%boundaries,fwi%k_max,           &
         fwi%epsilon,fwi%source_loc,fwi%receiver_loc,fwi%nb_frame,.true.)




    call free_acoustic_problem(fwi%forward)


  end subroutine one_fwi_step


end module m_fwi
