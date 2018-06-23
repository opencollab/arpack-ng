c
c   Flags for parallel execution
c   to be cleared on each new execution of parpack
c
c
      subroutine pcontext
      include 'pcontext.h'
      apps_first = .true.
      aitr_first = .true.
      end subroutine