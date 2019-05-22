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

typedef struct
{
  QPHIX_Real real;
  QPHIX_Real imag;
} QPHIX_Complex;

typedef struct
{ 
  QPHIX_Complex c[3];
} SU3_Vector;

typedef struct
{
  QPHIX_Complex e[3][3];
} SU3_Matrix;

typedef struct
{
  QPHIX_Complex m01, m02, m12;
  SU3_Real m00im, m11im, m22im;
  SU3_Real space;
} SU3_AntiHermitMat;

typedef struct 
{
  void *even;
  void *odd;
  QPHIX_evenodd_t parity;
} QPHIX_ColorMatrix;

typedef struct
{
  QPHIX_Complex t4[3][3][3][3];
} QPHIX_ColorTensor4;

typedef QPHIX_Real QSU3V[3][2][VECLEN];
typedef QPHIX_Real QSU3M[3][3][2][VECLEN];
typedef QPHIX_Complex qphix_su3_matrix[3][3];

#endif
