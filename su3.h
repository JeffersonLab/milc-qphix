#ifndef _SU3_H
#define _SU3_H

#if (QPHIX_PrecisionInt == 1)
typedef float SU3_Real;
typedef float QPHIX_Real;
#endif
#if (QPHIX_PrecisionInt == 2)
typedef double SU3_Real;
typedef double QPHIX_Real;
#endif

#define QPHIX_real(a) (a).real
#define QPHIX_imag(a) (a).imag
#define QPHIX_elem_M(a,ic,jc) ((a)[ic][jc])
#define Uelem(a,b) QPHIX_elem_M(*U,a,b)

#endif
