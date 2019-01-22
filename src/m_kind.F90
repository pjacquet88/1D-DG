module m_kind
  implicit none
  
  integer, parameter :: sp = selected_real_kind(6,37) ! single precision
  integer, parameter :: dp = selected_real_kind(15,307) ! double precision

  integer,parameter :: mp = selected_real_kind(KIND_REAL) ! my real precision

end module m_kind
