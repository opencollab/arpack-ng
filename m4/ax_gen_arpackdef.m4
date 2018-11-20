#
# SYNOPSIS
#
#   AX_GEN_ARPACKDEF()
#
# DESCRIPTION
#
#   This macro handles #define symbols for dedicated architecture.
#   In most cases, this is not necessary. For specific architecture
#   (like ILP64), some symbols must be changed accordingly.
#
#   Autoheader generates arpackdef.autotools.h (from arpackdef.autotools.h.in)
#   to handle #define symbols for C. The file arpackdef.autotools.h (C-style)
#   may not be used by fortran.
#
#   AX_GEN_ARPACKDEF generates arpackdef.h from arpackdef.autotools.h.
#   arpackdef.h can be used by fortran and/or C.
#

AC_DEFUN([AX_GEN_ARPACKDEF],
         [
           AC_MSG_NOTICE([generating arpackdef.h])
           echo "#ifndef __ARPACKDEF_H__"      >  arpackdef.h
           echo "#define __ARPACKDEF_H__"      >> arpackdef.h
           echo ""                             >> arpackdef.h
           grep ^.define arpackdef.autotools.h >> arpackdef.h
           echo ""                             >> arpackdef.h
           echo "#if INTERFACE64"              >> arpackdef.h
           echo "#define c_int  c_int64_t"     >> arpackdef.h
           echo "#define a_int    int64_t"     >> arpackdef.h
           echo "#define a_uint  uint64_t"     >> arpackdef.h
           echo "#else"                        >> arpackdef.h
           echo "#define a_int            int" >> arpackdef.h
           echo "#define a_uint  unsigned int" >> arpackdef.h
           echo "#endif"                       >> arpackdef.h
           echo ""                             >> arpackdef.h
           echo "#endif"                       >> arpackdef.h
         ])
