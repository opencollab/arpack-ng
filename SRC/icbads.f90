! icba : iso_c_binding for arpack

subroutine dsaupd_c(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
                    iparam, ipntr, workd, workl, lworkl, info)        &
                    bind(c, name="dsaupd_c")
  use :: iso_c_binding
  implicit none
  integer(kind=c_int),                         intent(inout) :: ido
  character(kind=c_char), dimension(1),        intent(in)    :: bmat
  integer(kind=c_int),    value,               intent(in)    :: n
  character(kind=c_char), dimension(2),        intent(in)    :: which
  integer(kind=c_int),    value,               intent(in)    :: nev
  real(kind=c_double),    value,               intent(in)    :: tol
  real(kind=c_double),    dimension(n),        intent(inout) :: resid
  integer(kind=c_int),    value,               intent(in)    :: ncv
  real(kind=c_double),    dimension(ldv, ncv), intent(out)   :: v
  integer(kind=c_int),    value,               intent(in)    :: ldv
  integer(kind=c_int),    dimension(11),       intent(inout) :: iparam
  integer(kind=c_int),    dimension(11),       intent(out)   :: ipntr
  real(kind=c_double),    dimension(3*n),      intent(out)   :: workd
  real(kind=c_double),    dimension(lworkl),   intent(out)   :: workl
  integer(kind=c_int),    value,               intent(in)    :: lworkl
  integer(kind=c_int),                         intent(inout) :: info
  call dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
              iparam, ipntr, workd, workl, lworkl, info)
end subroutine dsaupd_c

subroutine dseupd_c(rvec, howmny, select, d, z, ldz, sigma,      &
                    bmat, n, which, nev, tol, resid, ncv, v, ldv,&
                    iparam, ipntr, workd, workl, lworkl, info)   &
                    bind(c, name="dseupd_c")
  use :: iso_c_binding
  implicit none
  logical(kind=c_bool),   value,               intent(in)    :: rvec
  character(kind=c_char), dimension(1),        intent(in)    :: howmny
  logical(kind=c_bool),   dimension(ncv),      intent(in)    :: select
  real(kind=c_double),    dimension(nev),      intent(out)   :: d
  real(kind=c_double),    dimension(n, nev),   intent(out)   :: z
  integer(kind=c_int),    value,               intent(in)    :: ldz
  real(kind=c_double),    value,               intent(in)    :: sigma
  character(kind=c_char), dimension(1),        intent(in)    :: bmat
  integer(kind=c_int),    value,               intent(in)    :: n
  character(kind=c_char), dimension(2),        intent(in)    :: which
  integer(kind=c_int),    value,               intent(in)    :: nev
  real(kind=c_double),    value,               intent(in)    :: tol
  real(kind=c_double),    dimension(n),        intent(inout) :: resid
  integer(kind=c_int),    value,               intent(in)    :: ncv
  real(kind=c_double),    dimension(ldv, ncv), intent(out)   :: v
  integer(kind=c_int),    value,               intent(in)    :: ldv
  integer(kind=c_int),    dimension(11),       intent(inout) :: iparam
  integer(kind=c_int),    dimension(11),       intent(out)   :: ipntr
  real(kind=c_double),    dimension(3*n),      intent(out)   :: workd
  real(kind=c_double),    dimension(lworkl),   intent(out)   :: workl
  integer(kind=c_int),    value,               intent(in)    :: lworkl
  integer(kind=c_int),                         intent(inout) :: info
  call dseupd(rvec, howmny, select, d, z, ldz, sigma,      &
              bmat, n, which, nev, tol, resid, ncv, v, ldv,&
              iparam, ipntr, workd, workl, lworkl, info)
end subroutine dseupd_c
