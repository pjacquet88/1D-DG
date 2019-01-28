module m_animation
 use m_file_function,  only : str
  implicit none

  public  :: progress_bar,gif_creation
  private :: gif_creation_python,gif_creation_gnuplot

contains


  subroutine gif_creation(gnuplot,animation,frame_step,modeling_max_step,fwi_max_step)
    implicit none
    logical         ,intent(in) :: gnuplot
    character(len=*),intent(in) :: animation
    integer         ,intent(in) :: frame_step
    integer         ,intent(in) :: modeling_max_step
    integer         ,intent(in) :: fwi_max_step
    integer                     :: max_step

    if (animation.eq.'model_update') then
       max_step=fwi_max_step
    else
       max_step=modeling_max_step
    end if

    if (gnuplot) then
       call gif_creation_gnuplot(animation,frame_step,max_step)
    else
       call gif_creation_python(animation,frame_step,max_step)
    end if
  end subroutine gif_creation


  subroutine gif_creation_python(animation,frame_step,max_step)
    implicit none
    character(len=*),intent(in) :: animation
    integer         ,intent(in) :: frame_step
    integer         ,intent(in) :: max_step
    character(len=100)          :: command
    character(len=10)           :: prefix

    if (animation.eq.'no') then
       print*,'No animation has been made'
    else if (animation.eq.'model_update') then
       prefix='VP'
    elseif(animation.eq.'data_forward') then
       prefix='FP'
    elseif(animation.eq.'fwi_forward') then
       prefix='P'
    elseif(animation.eq.'fwi_backward') then
       prefix='QP'
    end if

    if (animation.ne.'no') then
       command='python3 ../animation_script/animation.py'//' '                  &
            //trim(prefix)//' '//trim(str(frame_step))//' '//trim(str(max_step))
       call system(command)
    end if

  end subroutine gif_creation_python


  subroutine gif_creation_gnuplot(animation,frame_step,max_step)
    implicit none
    character(len=*),intent(in) :: animation
    integer         ,intent(in) :: frame_step
    integer         ,intent(in) :: max_step

    if (animation.eq.'no') then
       print*,'No animation has been made'
    else if (animation.eq.'model_update') then

       open(unit=78,file='../animation_script/script.gnuplot',action='write')
       write(78,*)'load "../animation_script/trace1.gnuplot"'
       write(78,*)'n=',max_step
       write(78,*)'a=',frame_step
       write(78,*)'load "../animation_script/trace2.gnuplot"'
       close(78)

       open(unit=79,file='../animation_script/trace2.gnuplot',action='write')
       write(79,*)'dt=0.1/n'
       write(79,*)'i=0'
       write(79,*)'set yrange[0:2]'
       write(79,*)'load "../animation_script/animate.gnuplot"'
       close(79)

       open(unit=80,file='../animation_script/animate.gnuplot',action='write')
       write(80,*) "plot '../Files/VP'.i.'.dat' w l title 'Velocity model'.i"
       write(80,*) "i=i+a"
       write(80,*) "if (i<n) reread"
       close(80)

       call system("gnuplot ../animation_script/script.gnuplot")
       call system("eog animate.gif")

    elseif(animation.eq.'data_forward') then

       open(unit=78,file='../animation_script/script.gnuplot',action='write')
       write(78,*)'load "../animation_script/trace1.gnuplot"'
       write(78,*)'n=',max_step
       write(78,*)'a=',frame_step
       write(78,*)'load "../animation_script/trace2.gnuplot"'
       close(78)

       open(unit=79,file='../animation_script/trace2.gnuplot',action='write')
       write(79,*)'dt=0.1/n'
       write(79,*)'i=0'
       write(79,*)'set yrange[-0.5:0.5]'
       write(79,*)'load "../animation_script/animate.gnuplot"'
       close(79)

       open(unit=80,file='../animation_script/animate.gnuplot',action='write')
       write(80,*) "plot '../Files/FP'.i.'.dat' w l title 'Forward data pressure'.i"
       write(80,*) "i=i+a"
       write(80,*) "if (i<n) reread"
       close(80)

       call system("gnuplot ../animation_script/script.gnuplot")
       call system("eog animate.gif")

    elseif(animation.eq.'fwi_forward') then

       open(unit=78,file='../animation_script/script.gnuplot',action='write')
       write(78,*)'load "../animation_script/trace1.gnuplot"'
       write(78,*)'n=',max_step
       write(78,*)'a=',frame_step
       write(78,*)'load "../animation_script/trace2.gnuplot"'
       close(78)

       open(unit=79,file='../animation_script/trace2.gnuplot',action='write')
       write(79,*)'dt=0.1/n'
       write(79,*)'i=0'
       write(79,*)'set yrange[-0.5:0.5]'
       write(79,*)'load "../animation_script/animate.gnuplot"'
       close(79)

       open(unit=80,file='../animation_script/animate.gnuplot',action='write')
       write(80,*) "plot '../Files/P'.i.'.dat' w l title 'Forward fwi pressure'.i"
       write(80,*) "i=i+a"
       write(80,*) "if (i<n) reread"
       close(80)

       call system("gnuplot ../animation_script/script.gnuplot")
       call system("eog animate.gif")

    elseif(animation.eq.'fwi_backward') then

       open(unit=78,file='../animation_script/script.gnuplot',action='write')
       write(78,*)'load "../animation_script/trace1.gnuplot"'
       write(78,*)'n=',max_step
       write(78,*)'a=',frame_step
       write(78,*)'load "../animation_script/trace2.gnuplot"'
       close(78)

       open(unit=79,file='../animation_script/trace2.gnuplot',action='write')
       write(79,*)'dt=0.1/n'
       write(79,*)'i=0'
       write(79,*)'set yrange[-0.2:0.2]'
       write(79,*)'load "../animation_script/animate.gnuplot"'
       close(79)

       open(unit=80,file='../animation_script/animate.gnuplot',action='write')
       write(80,*) "plot '../Files/QP'.i.'.dat' w l title 'Backward fwi pressure'.i"
       write(80,*) "i=i+a"
       write(80,*) "if (i<n) reread"
       close(80)

       call system("gnuplot ../animation_script/script.gnuplot")
       call system("eog animate.gif")
    end if
  end subroutine gif_creation_gnuplot

subroutine progress_bar(iter,max_iter,o_bar_char,o_bar_length)
  implicit none
  integer,intent(in) :: iter     ! current iteration
  integer,intent(in) :: max_iter ! final iteration


  character,optional :: o_bar_char   ! optional for bar loading character
  integer  ,optional :: o_bar_length ! optional for progress bar length

  character,parameter :: default_bar_char='#'  ! default value for bar_char
  integer  ,parameter :: default_bar_length=50 ! default value for bar_length

  integer                            :: bar_length
  integer                            :: table_length
  character,dimension(:),allocatable :: table
  character                          :: bar_char
  character(len=30)                  :: fmt

  integer :: fraction,pourcentage,i

  if(present(o_bar_char)) then
     bar_char=o_bar_char
  else
     bar_char=default_bar_char
  end if

  if (present(o_bar_length)) then
     bar_length=o_bar_length
     table_length=bar_length+4
  else
     bar_length=default_bar_length
     table_length=bar_length+4
  end if

  allocate(table(table_length))

  fraction=floor(real(iter*bar_length/max_iter))
  pourcentage=floor(iter*100.0/max_iter)

  table(1)='['
  do i=2,fraction+1
     table(i)=bar_char
  end do
  do i=fraction+2,bar_length+1
     table(i)=' '
  end do
  table(bar_length+2)=']'
  table(bar_length+3)='-'
  table(bar_length+4)='['

  ! NB : write can't read and flush a string, we need a table of character
  write(fmt,"(a,i0,a)") "(1a1,(a,i0,a,i0,a,",table_length,"a,i0,a),$)"
  write(*,fmt)  char(13),'Iteration ',iter,'/',max_iter,' ',table,pourcentage,'%]'

  if (iter.eq.max_iter) then
     write(*,*)
  end if

end subroutine progress_bar

end module m_animation
