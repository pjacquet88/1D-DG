program main
  use m_kind
  use m_polynom
  implicit none

  integer  :: order,i,j,k
  real(mp) :: x,sum,error
  logical  :: bool

  ! ----- SET RANDOM SESSION
  integer                           :: values(1:8)
  integer,dimension(:), allocatable :: seed

  bool=.true.

  call date_and_time(values=values)
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)

  do order=1,10
     call init_polynom(order)
     do i=1,100
        sum=0.0_mp
        call random_number(x)

        do j=1,order+1
           sum=sum+eval_polynom_b(base_b(j,:),x)
        end do

        error=abs(1.0_mp-sum)

        if (mp.eq.sp) then
           bool=(bool.and.(error.lt.1e-6_mp))
        else if (mp.eq.dp) then
           bool=(bool.and.(error.lt.1e-15_mp))
        end if
     end do

     call my_return(bool)
     call free_polynom

  end do

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
