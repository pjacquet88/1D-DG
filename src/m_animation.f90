module m_animation
  implicit none

contains

  subroutine gif_creation(animation,nb_iter_fwi,data_time_step)
    implicit none
    character(len=*),intent(in) :: animation
    integer         ,intent(in) :: nb_iter_fwi
    integer         ,intent(in) :: data_time_step

    integer :: n_frame

    
    if (animation.eq.'no') then
       print*,'No animation has been made'
    else if (animation.eq.'model_update') then

       n_frame=nb_iter_fwi

       open(unit=78,file='../gnuplot_script/script.gnuplot',action='write')
       write(78,*)'load "../gnuplot_script/trace1.gnuplot"'
       write(78,*)'n=',nb_iter_fwi
       write(78,*)'a=',1
       write(78,*)'load "../gnuplot_script/trace2.gnuplot"'
       close(78)

       open(unit=79,file='../gnuplot_script/trace2.gnuplot',action='write')
       write(79,*)'dt=0.1/n'
       write(79,*)'i=0'
       write(79,*)'set yrange[0:2]'
       write(79,*)'load "../gnuplot_script/animate.gnuplot"'
       close(79)

       open(unit=80,file='../gnuplot_script/animate.gnuplot',action='write')
       write(80,*) "plot '../Files/VP'.i.'.dat' w l title 'Velocity model'.i"
       write(80,*) "i=i+a"
       write(80,*) "if (i<n) reread"
       close(80)

       call system("gnuplot ../gnuplot_script/script.gnuplot")
       call system("eog animate.gif")

    elseif(animation.eq.'data_forward') then

       n_frame=data_time_step

       open(unit=78,file='../gnuplot_script/script.gnuplot',action='write')
       write(78,*)'load "../gnuplot_script/trace1.gnuplot"'
       write(78,*)'n=',n_frame
       write(78,*)'a=',10
       write(78,*)'load "../gnuplot_script/trace2.gnuplot"'
       close(78)

       open(unit=79,file='../gnuplot_script/trace2.gnuplot',action='write')
       write(79,*)'dt=0.1/n'
       write(79,*)'i=0'
       write(79,*)'set yrange[-0.5:0.5]'
       write(79,*)'load "../gnuplot_script/animate.gnuplot"'
       close(79)

       open(unit=80,file='../gnuplot_script/animate.gnuplot',action='write')
       write(80,*) "plot '../Files/FP'.i.'.dat' w l title 'Forward data pressure'.i"
       write(80,*) "i=i+a"
       write(80,*) "if (i<n) reread"
       close(80)

       call system("gnuplot ../gnuplot_script/script.gnuplot")
       call system("eog animate.gif")

    elseif(animation.eq.'fwi_forward') then

       n_frame=data_time_step
       open(unit=78,file='../gnuplot_script/script.gnuplot',action='write')
       write(78,*)'load "../gnuplot_script/trace1.gnuplot"'
       write(78,*)'n=',n_frame
       write(78,*)'a=',10
       write(78,*)'load "../gnuplot_script/trace2.gnuplot"'
       close(78)

       open(unit=79,file='../gnuplot_script/trace2.gnuplot',action='write')
       write(79,*)'dt=0.1/n'
       write(79,*)'i=0'
       write(79,*)'set yrange[-0.5:0.5]'
       write(79,*)'load "../gnuplot_script/animate.gnuplot"'
       close(79)

       open(unit=80,file='../gnuplot_script/animate.gnuplot',action='write')
       write(80,*) "plot '../Files/P'.i.'.dat' w l title 'Forward fwi pressure'.i"
       write(80,*) "i=i+a"
       write(80,*) "if (i<n) reread"
       close(80)

       call system("gnuplot ../gnuplot_script/script.gnuplot")
       call system("eog animate.gif")

    elseif(animation.eq.'fwi_backward') then

       n_frame=data_time_step
       open(unit=78,file='../gnuplot_script/script.gnuplot',action='write')
       write(78,*)'load "../gnuplot_script/trace1.gnuplot"'
       write(78,*)'n=',n_frame
       write(78,*)'a=',10
       write(78,*)'load "../gnuplot_script/trace2.gnuplot"'
       close(78)

       open(unit=79,file='../gnuplot_script/trace2.gnuplot',action='write')
       write(79,*)'dt=0.1/n'
       write(79,*)'i=0'
       write(79,*)'set yrange[-0.2:0.2]'
       write(79,*)'load "../gnuplot_script/animate.gnuplot"'
       close(79)

       open(unit=80,file='../gnuplot_script/animate.gnuplot',action='write')
       write(80,*) "plot '../Files/QP'.i.'.dat' w l title 'Backward fwi pressure'.i"
       write(80,*) "i=i+a"
       write(80,*) "if (i<n) reread"
       close(80)

       call system("gnuplot ../gnuplot_script/script.gnuplot")
       call system("eog animate.gif")
    end if
end subroutine gif_creation


subroutine progress_bar(iter,max_iter)
  implicit none
  integer,intent(in) :: iter
  integer,intent(in) :: max_iter

  integer,parameter :: nb_char=50

  integer :: fraction,pourcentage,i
  character,dimension(54) :: table

  fraction=floor(real(iter*nb_char/max_iter))
  pourcentage=floor(iter*100.0/max_iter)

  table(1)='['
  do i=2,fraction+1
     table(i)='#'
  end do
  do i=fraction+2,nb_char+1
     table(i)=' '
  end do
  table(52)=']'
  table(53)='-'
  table(54)='['
  !  write(*,'(1a1,a,$)') char(13), table
  
!  write(*,'(1a1,(10a,a,a,54a,i0,a),$)')  char(13),'Iteration ','/',' ',table,pourcentage,']'
  write(*,'(1a1,(a,i0,a,i0,a,54a,i0,a),$)')  char(13),'Iteration ',iter,'/',max_iter,' ',table,pourcentage,'%]'

  if (iter.eq.max_iter) then
     write(*,*)
  end if
end subroutine progress_bar

end module m_animation
