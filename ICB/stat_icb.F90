! iso_c_binding : statistics

subroutine sstats_c() bind(c, name="sstats_c")
  use :: iso_c_binding
  implicit none
  call sstats()
end subroutine sstats_c

subroutine sstatn_c() bind(c, name="sstatn_c")
  use :: iso_c_binding
  implicit none
  call sstatn()
end subroutine sstatn_c

subroutine cstatn_c() bind(c, name="cstatn_c")
  use :: iso_c_binding
  implicit none
  call cstatn()
end subroutine cstatn_c

subroutine stat_c(  nopx_c,    nbx_c, nrorth_c, nitref_c, nrstrt_c,                    &
                  tsaupd_c, tsaup2_c, tsaitr_c, tseigt_c, tsgets_c, tsapps_c, tsconv_c,&
                  tnaupd_c, tnaup2_c, tnaitr_c, tneigh_c, tngets_c, tnapps_c, tnconv_c,&
                  tcaupd_c, tcaup2_c, tcaitr_c, tceigh_c, tcgets_c, tcapps_c, tcconv_c,&
                  tmvopx_c,  tmvbx_c, tgetv0_c, titref_c,  trvec_c)                    &
                  bind(c, name="stat_c")
  use :: iso_c_binding
  implicit none
#include "arpackicb.h"
  integer(kind=i_int), intent(out) ::   nopx_c,    nbx_c, nrorth_c, nitref_c, nrstrt_c
  real(kind=c_float),  intent(out) :: tsaupd_c, tsaup2_c, tsaitr_c, tseigt_c, tsgets_c, tsapps_c, tsconv_c,&
                                      tnaupd_c, tnaup2_c, tnaitr_c, tneigh_c, tngets_c, tnapps_c, tnconv_c,&
                                      tcaupd_c, tcaup2_c, tcaitr_c, tceigh_c, tcgets_c, tcapps_c, tcconv_c,&
                                      tmvopx_c,  tmvbx_c, tgetv0_c, titref_c,  trvec_c
  include 'stat.h'
    nopx_c =   nopx
     nbx_c =    nbx
  nrorth_c = nrorth
  nitref_c = nitref
  nrstrt_c = nrstrt
  tsaupd_c = tsaupd
  tsaup2_c = tsaup2
  tsaitr_c = tsaitr
  tseigt_c = tseigt
  tsgets_c = tsgets
  tsapps_c = tsapps
  tsconv_c = tsconv
  tnaupd_c = tnaupd
  tnaup2_c = tnaup2
  tnaitr_c = tnaitr
  tneigh_c = tneigh
  tngets_c = tngets
  tnapps_c = tnapps
  tnconv_c = tnconv
  tcaupd_c = tcaupd
  tcaup2_c = tcaup2
  tcaitr_c = tcaitr
  tceigh_c = tceigh
  tcgets_c = tcgets
  tcapps_c = tcapps
  tcconv_c = tcconv
  tmvopx_c = tmvopx
   tmvbx_c =  tmvbx
  tgetv0_c = tgetv0
  titref_c = titref
   trvec_c =  trvec
end subroutine stat_c
