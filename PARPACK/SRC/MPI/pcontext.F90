!
!   Flags for parallel execution
!   to be cleared on each new execution of parpack
!
!
      subroutine pcontext
#include "pcontextF90.h"
      apps_first = .true.
      aitr_first = .true.
      end subroutine
