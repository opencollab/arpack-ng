! iso_c_binding : seeds

subroutine iseed_c(seed1_c, seed2_c, seed3_c, seed4_c)&
                   bind(c, name="iseed_c")
  use :: iso_c_binding
  implicit none
#include "arpackdef.h"
  integer(kind=c_int), intent(in) :: seed1_c, seed2_c, seed3_c, seed4_c
  !
  call siseed(seed1_c, seed2_c, seed3_c, seed4_c)
end subroutine iseed_c
