! icbp : iso_c_binding for parpack

subroutine pznaupd_c(comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
                     iparam, ipntr, workd, workl, lworkl, rwork, info)       &
                     bind(c, name="pznaupd_c")
  use :: iso_c_binding
  implicit none
#include "arpackicb.h"
  integer(kind=i_int),              value,              intent(in)    :: comm
  integer(kind=i_int),                                  intent(inout) :: ido
  character(kind=c_char),                               intent(in)    :: bmat
  integer(kind=i_int),              value,              intent(in)    :: n
  character(kind=c_char),           dimension(2),       intent(in)    :: which
  integer(kind=i_int),              value,              intent(in)    :: nev
  real(kind=c_double),              value,              intent(in)    :: tol
  complex(kind=c_double_complex),   dimension(n),       intent(inout) :: resid
  integer(kind=i_int),              value,              intent(in)    :: ncv
  complex(kind=c_double_complex),   dimension(ldv, ncv),intent(out)   :: v
  integer(kind=i_int),              value,              intent(in)    :: ldv
  integer(kind=i_int),              dimension(11),      intent(inout) :: iparam
  integer(kind=i_int),              dimension(14),      intent(out)   :: ipntr
  complex(kind=c_double_complex),   dimension(3*n),     intent(out)   :: workd
  complex(kind=c_double_complex),   dimension(lworkl),  intent(out)   :: workl
  integer(kind=i_int),              value,              intent(in)    :: lworkl
  complex(kind=c_double_complex),   dimension(ncv),     intent(out)   :: rwork
  integer(kind=i_int),                                  intent(inout) :: info
  
  character(len=2):: w
  integer         :: i
  
  do i =1,2
      w(i:i) = which(i)
  end do
  
  call pznaupd(comm, ido, bmat, n, w, nev, tol, resid, ncv, v, ldv,&
               iparam, ipntr, workd, workl, lworkl, rwork, info)
end subroutine pznaupd_c

subroutine pzneupd_c(comm, rvec, howmny, select, d, z, ldz, sigma, workev,&
                     bmat, n, which, nev, tol, resid, ncv, v, ldv,        &
                     iparam, ipntr, workd, workl, lworkl, rwork, info)    &
                     bind(c, name="pzneupd_c")
  use :: iso_c_binding
  implicit none
#include "arpackicb.h"
  integer(kind=i_int),              value,              intent(in)    :: comm
  integer(kind=i_int),              value,              intent(in)    :: rvec
  character(kind=c_char),                               intent(in)    :: howmny
  integer(kind=i_int),              dimension(ncv),     intent(in)    :: select
  complex(kind=c_double_complex),   dimension(nev),     intent(out)   :: d
  complex(kind=c_double_complex),   dimension(n, nev),  intent(out)   :: z
  integer(kind=i_int),              value,              intent(in)    :: ldz
  complex(kind=c_double_complex),   value,              intent(in)    :: sigma
  complex(kind=c_double_complex),   dimension(2*ncv),   intent(out)   :: workev
  character(kind=c_char),                               intent(in)    :: bmat
  integer(kind=i_int),              value,              intent(in)    :: n
  character(kind=c_char),           dimension(2),       intent(in)    :: which
  integer(kind=i_int),              value,              intent(in)    :: nev
  real(kind=c_double),              value,              intent(in)    :: tol
  complex(kind=c_double_complex),   dimension(n),       intent(inout) :: resid
  integer(kind=i_int),              value,              intent(in)    :: ncv
  complex(kind=c_double_complex),   dimension(ldv, ncv),intent(out)   :: v
  integer(kind=i_int),              value,              intent(in)    :: ldv
  integer(kind=i_int),              dimension(11),      intent(inout) :: iparam
  integer(kind=i_int),              dimension(14),      intent(out)   :: ipntr
  complex(kind=c_double_complex),   dimension(3*n),     intent(out)   :: workd
  complex(kind=c_double_complex),   dimension(lworkl),  intent(out)   :: workl
  integer(kind=i_int),              value,              intent(in)    :: lworkl
  complex(kind=c_double_complex),   dimension(ncv),     intent(out)   :: rwork
  integer(kind=i_int),                                  intent(inout) :: info

  ! convert parameters if needed.

  logical :: rv
  logical, dimension(ncv) :: slt
  integer :: idx
  character(len=2):: w
  integer         :: i
  
  rv = .false.
  if (rvec .ne. 0) rv = .true.

  slt = .false.
  do idx=1, ncv
    if (select(idx) .ne. 0) slt(idx) = .true.
  enddo
  
  do i =1,2
      w(i:i) = which(i)
  end do

  ! call arpack.

  call pzneupd(comm, rv, howmny, slt, d, z, ldz, sigma, workev,&
               bmat, n, w, nev, tol, resid, ncv, v, ldv,   &
               iparam, ipntr, workd, workl, lworkl, rwork, info)
end subroutine pzneupd_c
