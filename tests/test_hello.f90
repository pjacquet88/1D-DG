program main
  implicit none

  print*,'Hello World'

  call my_return(.true.)

contains

  subroutine my_return(bool)
    logical,intent(in) :: bool

    if (bool) then
       print*,'Succeed'
    else
       ERROR STOP "TEST FAILED"
    end if
  end subroutine my_return
end program main
