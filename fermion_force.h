#ifndef FERMION_FORCE_H
#define FERMION_FORCE_H

extern "C"
{
#include "qphix.h"
#include "su3.h"

#include "veclength.h"

/* Directions */
#define NDIRS          8           /* number of directions */
#define NODIR          -1          /* not a direction */
#define FORWARDS       1
#define BACKWARDS      (-1)        /* BACKWARDS = -FORWARDS */
#define OPP_DIR(dir)   (7-(dir))   /* Opposite direction */
#define XUP            0
#define YUP            1
#define ZUP            2
#define TUP            3
#define TDOWN          4
#define ZDOWN          5
#define YDOWN          6
#define XDOWN          7
/* defines for 3rd nearest neighbor (NAIK) stuff */
#define X3UP           8
#define Y3UP           9
#define Z3UP           10
#define T3UP           11
#define T3DOWN         12
#define Z3DOWN         13
#define Y3DOWN         14
#define X3DOWN         15
#define OPP_3_DIR(dir) (23-(dir))
#define DIR3(dir)      ((dir)+8)

#define QPHIX_MAX_NAIK 20
#define QPHIX_PI 3.14159265358979323846
#define QPHIX_EPS_NAIK_ZERO {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}

#if (QPHIX_PrecisionInt == 1)
#define QPHIX_U3_UNIT_ANALYTIC_EPS 1.0e-6
#define QPHIX_SVD3x3_PREC 5e-16
#endif
#if (QPHIX_PrecisionInt == 2)
#define QPHIX_U3_UNIT_ANALYTIC_EPS 1.0e-14
#define QPHIX_SVD3x3_PREC 5e-7
#endif

//typedef enum
//{
//  QPHIX_UNITARIZE_U3 = 0,
//  QPHIX_UNITARIZE_SU3 = 1,
//} QPHIX_hisq_unitarize_group_t;
//
//typedef enum
//{
//  QPHIX_UNITARIZE_NONE = 0,
//  QPHIX_UNITARIZE_RATIONAL = 3,
//  QPHIX_UNITARIZE_ANALYTIC = 5,
//} QPHIX_hisq_unitarize_method_t;

typedef struct
{
  int x;
  int y;
  int z;
  int t;
} site;

//void QPHIX_hisq_force_multi(QPHIX_info_t *info, QPHIX_FermionLinksHisq *flh,
//                            QPHIX_ColorMatrix *force[],
//                            QPHIX_hisq_coeffs_t *hisq_coeff,
//                            QPHIX_Real *residues, QPHIX_ColorVector *x[],
//                            int *n_orders_naik);

}

#endif
