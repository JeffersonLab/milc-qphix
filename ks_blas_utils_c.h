#ifndef _KS_BLAS_UTILS_C_H_
#define _KS_BLAS_UTILS_C_H_

#include "ks_globals.h"
#include "qcd_data_types.h" 

/* Non-optimized version of BLAS routines used in MILC's CG solver. Once these 
 * versions meet the correctness criteria, have to look at incorporating the 
 * versions already available in QPhiX.
 */

double
magsq_su3vec (const KS *spinor);
 
void 
scalar_mult_add_su3_vector (const KS *sp1, const KS *sp2, fptype s, KS *sp3);

void
scalar_mult_su3_vector (const KS *in, double scalar, KS *out);

double
su3_rdot (const KS *sp1, const KS *sp2);

void 
add_su3_vector (const KS *sp1, const KS *sp2, KS *sp3);

void
sub_su3_vector (const KS *sp1, const KS *sp2, KS *sp3);
 
void
g_doublesum (double *d);

fptype 
relative_residue (KS *p, KS *q, int parity);

/******************************** NEW BLAS API ********************************/

inline void 
scalar_mult_add_su3_vector2 (fptype in1, fptype in2, fptype s, fptype &out3)
{
    out3 = in1 + s*in2;
}

inline void
scalar_mult_su3_vector2 (fptype in, double scalar, fptype &out)
{
    out = in * scalar;
}

/* Add two KS* spinors */
inline void
add_su3_vector2 (fptype in1, fptype in2, fptype &out) 
{
    out = in1 + in2;
}

/* Subtract two KS* spinors */
inline void
sub_su3_vector2 (fptype in1, fptype in2, fptype &out)  
{
    out = in1 - in2;
}

/******************************************************************************/


#endif
