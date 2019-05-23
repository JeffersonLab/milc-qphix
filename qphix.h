#ifndef _QPHIX_H
#define _QPHIX_H

/* This file is for external access */

  // allow the user to specify QPHIX_Precision and/or QPHIX_PrecisionInt
  // default to single precision.  This macro is used to convert generic
  // names for types and procedures to precision-specific names
  // This macro can be defined for the entire compilation or for
  // individual compilation units.  It must be defined before including 
  // this file.
  // This logic is the same as the QOP logic.  It is a bit complicated
  // because it allows two waysto define the macro 
  // and it sets a default of single precision.

#ifndef QPHIX_Precision
#  ifndef QPHIX_PrecisionInt
//# warning "QPHIX_PrecisionInt NOT Defined"
#    define QPHIX_PrecisionInt 1
#else
//# warning "QPHIX_PrecisionInt Defined"
#  endif
#  if QPHIX_PrecisionInt == 1
#    define QPHIX_Precision 'F'
#    define QPHIX_PrecisionLetter F
#  elif QPHIX_PrecisionInt == 2
#    define QPHIX_Precision 'D'
#    define QPHIX_PrecisionLetter D
#  else
#    error "bad QPHIX_PrecisionInt"
#  endif
#else
#  ifndef QPHIX_PrecisionInt
#    if QPHIX_Precision == 'F'
#      define QPHIX_PrecisionInt 1
#      define QPHIX_PrecisionLetter F
#    elif QPHIX_Precision == 'D'
#      define QPHIX_PrecisionInt 2
#      define QPHIX_PrecisionLetter D
#    else
#      error "bad QPHIX_Precision"
#    endif
#  else
#    if QPHIX_Precision == 'F'
#      if QPHIX_PrecisionInt != 1
#        error "inconsistent QPHIX_Precision='F' and QPHIX_PrecisionInt"
#      endif
#      define QPHIX_PrecisionLetter F
#    elif QPHIX_Precision == 'D'
#      if QPHIX_PrecisionInt != 2
#        error "inconsistent QPHIX_Precision='D' and QPHIX_PrecisionInt"
#      endif
#      define QPHIX_PrecisionLetter D
#    else
#      error "bad QPHIX_Precision"
#    endif
#  endif
#endif

#include <qphix_int.h>
#include <qphix_f3.h>
#include <qphix_d3.h>
//#include <qphix_misc.h>
#endif
