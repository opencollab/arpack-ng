! icba : iso_c_binding for arpack

subroutine cnaupd_c(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
                    iparam, ipntr, workd, workl, lworkl, rwork, info) &
                    bind(c, name="cnaupd_c")
  use :: iso_c_binding
  implicit none
#include "arpackdef.h"
  integer(kind=c_int),                              intent(inout) :: ido
  character(kind=c_char),       dimension(1),       intent(in)    :: bmat
  integer(kind=c_int),          value,              intent(in)    :: n
  character(kind=c_char),       dimension(2),       intent(in)    :: which
  integer(kind=c_int),          value,              intent(in)    :: nev
  real(kind=c_float),           value,              intent(in)    :: tol
  complex(kind=c_float_complex),dimension(n),       intent(inout) :: resid
  integer(kind=c_int),          value,              intent(in)    :: ncv
  complex(kind=c_float_complex),dimension(ldv, ncv),intent(out)   :: v
  integer(kind=c_int),          value,              intent(in)    :: ldv
  integer(kind=c_int),          dimension(11),      intent(inout) :: iparam
  integer(kind=c_int),          dimension(11),      intent(out)   :: ipntr
  complex(kind=c_float_complex),dimension(3*n),     intent(out)   :: workd
  complex(kind=c_float_complex),dimension(lworkl),  intent(out)   :: workl
  integer(kind=c_int),          value,              intent(in)    :: lworkl
  real(kind=c_float),           dimension(ncv),     intent(out)   :: rwork
  integer(kind=c_int),                              intent(inout) :: info
  call cnaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
              iparam, ipntr, workd, workl, lworkl, rwork, info)
end subroutine cnaupd_c

subroutine cneupd_c(rvec, howmny, select, d, z, ldz, sigma, workev,  &
                    bmat, n, which, nev, tol, resid, ncv, v, ldv,    &
                    iparam, ipntr, workd, workl, lworkl, rwork, info)&
                    bind(c, name="cneupd_c")
  use :: iso_c_binding
  implicit none
#include "arpackdef.h"
  integer(kind=c_int),          value,              intent(in)    :: rvec
  character(kind=c_char),       dimension(1),       intent(in)    :: howmny
  integer(kind=c_int),          dimension(ncv),     intent(in)    :: select
  complex(kind=c_float_complex),dimension(nev),     intent(out)   :: d
  complex(kind=c_float_complex),dimension(n, nev),  intent(out)   :: z
  integer(kind=c_int),          value,              intent(in)    :: ldz
  complex(kind=c_float_complex),value,              intent(in)    :: sigma
  complex(kind=c_float_complex),dimension(2*ncv),   intent(out)   :: workev
  character(kind=c_char),       dimension(1),       intent(in)    :: bmat
  integer(kind=c_int),          value,              intent(in)    :: n
  character(kind=c_char),       dimension(2),       intent(in)    :: which
  integer(kind=c_int),          value,              intent(in)    :: nev
  real(kind=c_float),           value,              intent(in)    :: tol
  complex(kind=c_float_complex),dimension(n),       intent(inout) :: resid
  integer(kind=c_int),          value,              intent(in)    :: ncv
  complex(kind=c_float_complex),dimension(ldv, ncv),intent(out)   :: v
  integer(kind=c_int),          value,              intent(in)    :: ldv
  integer(kind=c_int),          dimension(11),      intent(inout) :: iparam
  integer(kind=c_int),          dimension(11),      intent(out)   :: ipntr
  complex(kind=c_float_complex),dimension(3*n),     intent(out)   :: workd
  complex(kind=c_float_complex),dimension(lworkl),  intent(out)   :: workl
  integer(kind=c_int),          value,              intent(in)    :: lworkl
  real(kind=c_float),           dimension(ncv),     intent(out)   :: rwork
  integer(kind=c_int),                              intent(inout) :: info

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

  call cneupd(rv, howmny, slt, d, z, ldz, sigma, workev,     &
              bmat, n, which, nev, tol, resid, ncv, v, ldv,  &
              iparam, ipntr, workd, workl, lworkl, rwork, info)
end subroutine cneupd_c
