! icbp : iso_c_binding for parpack

subroutine psnaupd_c(comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
                     iparam, ipntr, workd, workl, lworkl, info)              &
                     bind(c, name="psnaupd_c")
  use :: iso_c_binding
  implicit none
#include "arpackdef.h"
  integer(kind=c_int),    value,               intent(in)    :: comm
  integer(kind=c_int),                         intent(inout) :: ido
  character(kind=c_char), dimension(1),        intent(in)    :: bmat
  integer(kind=c_int),    value,               intent(in)    :: n
  character(kind=c_char), dimension(2),        intent(in)    :: which
  integer(kind=c_int),    value,               intent(in)    :: nev
  real(kind=c_float),     value,               intent(in)    :: tol
  real(kind=c_float),     dimension(n),        intent(inout) :: resid
  integer(kind=c_int),    value,               intent(in)    :: ncv
  real(kind=c_float),     dimension(ldv, ncv), intent(out)   :: v
  integer(kind=c_int),    value,               intent(in)    :: ldv
  integer(kind=c_int),    dimension(11),       intent(inout) :: iparam
  integer(kind=c_int),    dimension(11),       intent(out)   :: ipntr
  real(kind=c_float),     dimension(3*n),      intent(out)   :: workd
  real(kind=c_float),     dimension(lworkl),   intent(out)   :: workl
  integer(kind=c_int),    value,               intent(in)    :: lworkl
  integer(kind=c_int),                         intent(inout) :: info
  call psnaupd(comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
               iparam, ipntr, workd, workl, lworkl, info)
end subroutine psnaupd_c

subroutine psneupd_c(comm, rvec, howmny, select,                  &
                     dr, di, z, ldz, sigmar, sigmai, workev,      &
                     bmat, n, which, nev, tol, resid, ncv, v, ldv,&
                     iparam, ipntr, workd, workl, lworkl, info)   &
                     bind(c, name="psneupd_c")
  use :: iso_c_binding
  implicit none
#include "arpackdef.h"
  integer(kind=c_int),    value,               intent(in)    :: comm
  integer(kind=c_int),    value,               intent(in)    :: rvec
  character(kind=c_char), dimension(1),        intent(in)    :: howmny
  integer(kind=c_int),    dimension(ncv),      intent(in)    :: select
  real(kind=c_float),     dimension(nev+1),    intent(out)   :: dr
  real(kind=c_float),     dimension(nev+1),    intent(out)   :: di
  real(kind=c_float),     dimension(n, nev+1), intent(out)   :: z
  integer(kind=c_int),    value,               intent(in)    :: ldz
  real(kind=c_float),     value,               intent(in)    :: sigmar
  real(kind=c_float),     value,               intent(in)    :: sigmai
  real(kind=c_float),     dimension(3*ncv),    intent(out)   :: workev
  character(kind=c_char), dimension(1),        intent(in)    :: bmat
  integer(kind=c_int),    value,               intent(in)    :: n
  character(kind=c_char), dimension(2),        intent(in)    :: which
  integer(kind=c_int),    value,               intent(in)    :: nev
  real(kind=c_float),     value,               intent(in)    :: tol
  real(kind=c_float),     dimension(n),        intent(inout) :: resid
  integer(kind=c_int),    value,               intent(in)    :: ncv
  real(kind=c_float),     dimension(ldv, ncv), intent(out)   :: v
  integer(kind=c_int),    value,               intent(in)    :: ldv
  integer(kind=c_int),    dimension(11),       intent(inout) :: iparam
  integer(kind=c_int),    dimension(11),       intent(out)   :: ipntr
  real(kind=c_float),     dimension(3*n),      intent(out)   :: workd
  real(kind=c_float),     dimension(lworkl),   intent(out)   :: workl
  integer(kind=c_int),    value,               intent(in)    :: lworkl
  integer(kind=c_int),                         intent(inout) :: info

  ! convert parameters if needed.

  logical :: rv
  logical, dimension(ncv) :: slt
  integer :: idx

  rv = .false.
  if (rvec .ne. 0) rv = .true.

  slt = .false.
  do idx=1, ncv
    if (select(idx) .ne. 0) slt(idx) = .true.
  enddo

  ! call arpack.

  call psneupd(comm, rv, howmny, slt,                       &
               dr, di, z, ldz, sigmar, sigmai, workev,      &
               bmat, n, which, nev, tol, resid, ncv, v, ldv,&
               iparam, ipntr, workd, workl, lworkl, info)
end subroutine psneupd_c
