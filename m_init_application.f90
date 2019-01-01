module m_init_application

  implicit none

    !*************** Problem Parameters *******************************************
  integer           :: nb_elem        ! Nb of elements (all same length)
  integer           :: order,DoF      ! Polynoms order
  real              :: total_length   ! domain length
  real              :: final_time     ! final time
  character(len=20) :: time_scheme    ! change the time scheme
  real              :: alpha          ! Penalisation value
  character(len=20) :: signal         ! initial values (flat = 0)
  character(len=20) :: boundaries     ! Boundary Conditions
  logical           :: bernstein      ! If F-> Lagrange Elements
  integer           :: k_max          ! iterr max for power method algo.
  real              :: epsilon        ! precision for power method algo.
  integer           :: source_loc     ! location of the source (elemts)
  integer           :: receiver_loc   ! location of the receiver(elemts)
  character(len=20) :: strategy
  logical           :: adjoint_test
  character(len=20) :: scalar_product
  integer           :: nb_iter_fwi

  real,dimension(:),allocatable :: velocity_data,velocity_ini
  real,dimension(:),allocatable :: density_data,density_ini
  integer                       :: size_velocity_data,size_velocity_ini
  integer                       :: size_density_data,size_density_ini

  character(len=20) :: animation


contains

  subroutine setup_parameters()
    implicit none
    integer           :: n_arg
    character(len=20) :: param_file_name
    n_arg=command_argument_count()

    if (n_arg.lt.1) then
       print*,''//achar(27)//'[31m NO PARAMETER FILES GIVEN '//achar(27)//'[0m'
       print*,''//achar(27)//'[31m DEFAULT VALUES WILL BE TAKEN '//achar(27)//'[0m'
       call default_value_init()
       call print_selected_parameters()
    else if (n_arg.eq.1) then
       call get_command_argument(n_arg,param_file_name)
       print*,''//achar(27)//'[32m PARAMETER FILE GIVEN : '//achar(27)//'[0m',param_file_name
       call read_param(param_file_name)
       call print_selected_parameters()
    else
       print*,'NUMBER OF ARGUMENT NOT VALID'
    end if

  end subroutine setup_parameters


  subroutine read_param(param_file_name)
    implicit none
    character(len=20),intent(in) :: param_file_name
    character(len=20) :: dummy

    open(unit=1,file=param_file_name)
    read(1,*) dummy                              ! 1
    read(1,*) dummy                              ! 2
    read(1,*) dummy                              ! 3
    read(1,*) dummy                              ! 4
    read(1,*) dummy                              ! 5
    read(1,*) total_length                       ! 6
    read(1,*) dummy                              ! 7
    read(1,*) final_time                         ! 8
    read(1,*) dummy                              ! 9
    read(1,*) dummy                              ! 10
    read(1,*) dummy                              ! 11
    read(1,*) nb_elem                            ! 12
    read(1,*) dummy                              ! 13
    read(1,*) order                              ! 14
    DoF=order+1
    read(1,*) dummy                              ! 15
    read(1,*) bernstein                          ! 16
    read(1,*) dummy                              ! 17
    read(1,*) alpha                              ! 18
    read(1,*) dummy                              ! 19
    read(1,*) time_scheme                        ! 20
    read(1,*) dummy                              ! 21
    read(1,*) dummy                              ! 22
    read(1,*) dummy                              ! 23
    read(1,*) signal                             ! 24
    read(1,*) dummy                              ! 25
    read(1,*) boundaries                         ! 26
    read(1,*) dummy                              ! 27
    read(1,*) dummy                              ! 28
    read(1,*) dummy                              ! 29
    read(1,*) nb_iter_fwi                        ! 30
    read(1,*) dummy                              ! 31
    read(1,*) strategy                           ! 32
    read(1,*) dummy                              ! 33
    read(1,*) scalar_product                     ! 34
    read(1,*) dummy                              ! 35
    read(1,*) adjoint_test                       ! 36
    read(1,*) dummy                              ! 37
    read(1,*) dummy                              ! 38
    read(1,*) dummy                              ! 39
    read(1,*) source_loc                         ! 40
    read(1,*) dummy                              ! 41
    read(1,*) receiver_loc                       ! 42
    read(1,*) dummy                              ! 43
    read(1,*) dummy                              ! 44
    read(1,*) dummy                              ! 45
    read(1,*) size_velocity_data                 ! 46
    allocate(velocity_data(size_velocity_data))
    read(1,*) dummy                              ! 47
    read(1,*) velocity_data                      ! 48
    read(1,*) dummy                              ! 49
    read(1,*) size_velocity_ini                  ! 50
    allocate(velocity_ini(size_velocity_ini))
    read(1,*) dummy                              ! 51
    read(1,*) velocity_ini                       ! 52
    read(1,*) dummy                              ! 53 
    read(1,*) dummy                              ! 54
    read(1,*) size_density_data                  ! 55
    allocate(density_data(size_density_data))
    read(1,*) dummy                              ! 56
    read(1,*) density_data                       ! 57
    read(1,*) dummy                              ! 58
    read(1,*) size_density_ini                   ! 59
    allocate(density_ini(size_density_data))
    read(1,*) dummy                              ! 60
    read(1,*) density_ini                        ! 61
    read(1,*) dummy                              ! 62
    read(1,*) dummy                              ! 63
    read(1,*) dummy                              ! 64
    read(1,*) animation                          ! 65

    k_max=1e3
    epsilon=1e-5

  end subroutine read_param

  subroutine default_value_init()
    implicit none

    total_length=1.0
    final_time=3.0
    nb_elem=100
    order=2
    DoF=order+1
    bernstein=.true.
    alpha=1.0
    time_scheme='RK4'
    signal='flat'
    boundaries='ABC'
    nb_iter_fwi=25
    strategy='ATD'
    scalar_product='canonical'
    adjoint_test=.false.
    source_loc=1
    receiver_loc=30
    size_velocity_data=1
    allocate(velocity_data(size_velocity_data))
    velocity_data=1.0
    size_velocity_ini=1
    allocate(velocity_ini(size_velocity_ini))
    velocity_ini=1.0
    size_density_data=1
    allocate(density_data(size_density_data))
    density_data=1.0
    size_density_ini=1
    allocate(density_ini(size_density_ini))
    density_ini=1.0
    animation='data_forward'
    k_max=1e3
    epsilon=1e-5

  end subroutine default_value_init


  subroutine print_selected_parameters()
    implicit none

    print*,''//achar(27)//'[96m####################### PARAMETERS SELECTED ############################ '//achar(27)//'[0m'
    print*,''
    print*,''//achar(27)//'[96m####################### PROBLEM SIZE ################################### '//achar(27)//'[0m'
    print*,''//achar(27)//'[92mtotal_length       ='//achar(27)//'[0m',total_length
    print*,''//achar(27)//'[92mfinal_time         ='//achar(27)//'[0m',final_time
    print*,''//achar(27)//'[96m####################### PROBLEM DISCRETIZATION ######################### '//achar(27)//'[0m'
    print*,''//achar(27)//'[92mnb_elem            ='//achar(27)//'[0m',nb_elem
    print*,''//achar(27)//'[92morder              ='//achar(27)//'[0m',order
    print*,''//achar(27)//'[92mDoF                ='//achar(27)//'[0m',DoF
    print*,''//achar(27)//'[92mbernstein          ='//achar(27)//'[0m',bernstein
    print*,''//achar(27)//'[92malpha              ='//achar(27)//'[0m',alpha
    print*,''//achar(27)//'[92mtime_scheme        ='//achar(27)//'[0m',time_scheme
    print*,''//achar(27)//'[96m####################### DIRECT PROBLEM SPECIFICITIES ###################'//achar(27)//'[0m'
    print*,''//achar(27)//'[92msignal             ='//achar(27)//'[0m',signal
    print*,''//achar(27)//'[92mboundaries         ='//achar(27)//'[0m',boundaries
    print*,''//achar(27)//'[96m###################### FWI SPECIFICITIES ###############################'//achar(27)//'[0m'
    print*,''//achar(27)//'[92mnb_iter_fwi        ='//achar(27)//'[0m',nb_iter_fwi
    print*,''//achar(27)//'[92mstrategy           ='//achar(27)//'[0m',strategy
    print*,''//achar(27)//'[92mscalar_product     ='//achar(27)//'[0m',scalar_product
    print*,''//achar(27)//'[92madjoint_test       ='//achar(27)//'[0m',adjoint_test
    print*,''//achar(27)//'[96m#################### SOURCE AND RECEIVER LOCATION ######################'//achar(27)//'[0m'
    print*,''//achar(27)//'[92msource_loc         ='//achar(27)//'[0m',source_loc
    print*,''//achar(27)//'[92mreceiver_loc       ='//achar(27)//'[0m',receiver_loc
    print*,''//achar(27)//'[96m##################### VELOCITY MODEL ###################################'//achar(27)//'[0m'
    print*,''//achar(27)//'[92msize_velocity_data ='//achar(27)//'[0m',size_velocity_data
    print*,''//achar(27)//'[92mvelocity_data      ='//achar(27)//'[0m',velocity_data
    print*,''//achar(27)//'[92msize_velocity_ini  ='//achar(27)//'[0m',size_velocity_ini
    print*,''//achar(27)//'[92mvelocity_ini       ='//achar(27)//'[0m',velocity_ini
    print*,''//achar(27)//'[96m##################### DENSITY MODEL ####################################'//achar(27)//'[0m'
    print*,''//achar(27)//'[92msize_density_data  ='//achar(27)//'[0m',size_density_data
    print*,''//achar(27)//'[92mdensity_data       ='//achar(27)//'[0m',density_data
    print*,''//achar(27)//'[92msize_density_ini   ='//achar(27)//'[0m',size_density_ini
    print*,''//achar(27)//'[92mdensity_ini        ='//achar(27)//'[0m',density_ini
    print*,''//achar(27)//'[96m##################### ANIMATION ########################################'//achar(27)//'[0m'
    print*,''//achar(27)//'[92manimation           ='//achar(27)//'[0m',animation
    print*,''
    print*,''
  end subroutine print_selected_parameters

end module m_init_application
