! iso_c_binding : debug

subroutine debug_c(logfil_c, ndigit_c, mgetv0_c,                                        &
                   msaupd_c, msaup2_c, msaitr_c, mseigt_c, msapps_c, msgets_c, mseupd_c,&
                   mnaupd_c, mnaup2_c, mnaitr_c, mneigh_c, mnapps_c, mngets_c, mneupd_c,&
                   mcaupd_c, mcaup2_c, mcaitr_c, mceigh_c, mcapps_c, mcgets_c, mceupd_c)&
                   bind(c, name="debug_c")
  use :: iso_c_binding
  implicit none
#include "arpackicb.h"
  integer(kind=i_int), value, intent(in) :: logfil_c, ndigit_c, mgetv0_c
  integer(kind=i_int), value, intent(in) :: msaupd_c, msaup2_c, msaitr_c, mseigt_c, msapps_c, msgets_c, mseupd_c
  integer(kind=i_int), value, intent(in) :: mnaupd_c, mnaup2_c, mnaitr_c, mneigh_c, mnapps_c, mngets_c, mneupd_c
  integer(kind=i_int), value, intent(in) :: mcaupd_c, mcaup2_c, mcaitr_c, mceigh_c, mcapps_c, mcgets_c, mceupd_c
  include 'debug.h'
  logfil = logfil_c
  ndigit = ndigit_c
  mgetv0 = mgetv0_c
  msaupd = msaupd_c
  msaup2 = msaup2_c
  msaitr = msaitr_c
  mseigt = mseigt_c
  msapps = msapps_c
  msgets = msgets_c
  mseupd = mseupd_c
  mnaupd = mnaupd_c
  mnaup2 = mnaup2_c
  mnaitr = mnaitr_c
  mneigh = mneigh_c
  mnapps = mnapps_c
  mngets = mngets_c
  mneupd = mneupd_c
  mcaupd = mcaupd_c
  mcaup2 = mcaup2_c
  mcaitr = mcaitr_c
  mceigh = mceigh_c
  mcapps = mcapps_c
  mcgets = mcgets_c
  mceupd = mceupd_c
end subroutine debug_c
