#ifndef QPHIX_SU3_ALGEBRA_H
#define QPHIX_SU3_ALGEBRA_H

extern "C"
{
#include "fermion_force.h"

// Scalar version of SU3 algebra primitives
void QSU3_M_eq_Ma(qphix_su3_matrix *restrict r, qphix_su3_matrix *restrict a);

void QSU3_M_eq_zero ( qphix_su3_matrix *restrict r);

void QSU3_M_eq_Ma_times_M ( qphix_su3_matrix *restrict r, qphix_su3_matrix *restrict a, qphix_su3_matrix *restrict b);

void QSU3_M_eq_M_times_M ( qphix_su3_matrix *restrict r, qphix_su3_matrix *restrict a, qphix_su3_matrix *restrict b);

void QSU3_M_eq_M_times_Ma ( qphix_su3_matrix *restrict r, qphix_su3_matrix *restrict a, qphix_su3_matrix *restrict b);

void QSU3_R_eq_re_trace_M ( QPHIX_Real *restrict r, qphix_su3_matrix *restrict a);

void QSU3_M_eq_r_times_M ( qphix_su3_matrix *restrict r, QPHIX_Real *restrict a, qphix_su3_matrix *restrict b);

void QSU3_M_eq_r_times_M_plus_M ( qphix_su3_matrix *restrict r, QPHIX_Real *restrict a, qphix_su3_matrix *restrict b, qphix_su3_matrix *restrict c);


// Vectorized, multi-threaded, single processor version of SU3 alrebra primitves
void QPHIX_M_eq_zero(QPHIX_ColorMatrix *);

void QPHIX_M_eq_M(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *);

void QPHIX_M_peq_M(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *);

void QPHIX_M_eq_r_times_M(QPHIX_ColorMatrix *, QPHIX_Real, QPHIX_ColorMatrix *);

void QPHIX_M_peq_r_times_M(QPHIX_ColorMatrix *, QPHIX_Real, QPHIX_ColorMatrix *);

void QPHIX_M_peq_M_times_M(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *, QPHIX_ColorMatrix *);

void QPHIX_M_eq_M_times_M(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2);

void QPHIX_M_peq_M_times_Ma(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *, QPHIX_ColorMatrix *);

void QPHIX_M_eq_M_times_Ma(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2);

void QPHIX_M_peq_Ma(QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2);

void QPHIX_M_eq_Ma(QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2);

void QPHIX_M_eq_Ma_times_M(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2);

void QPHIX_M_eq_r_times_M_plus_M(QPHIX_ColorMatrix *M3, QPHIX_Real r, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2);

void QPHIX_M_eq_antiherm_M( QPHIX_ColorMatrix *m1, QPHIX_ColorMatrix *m2 );

void QPHIX_M_eqm_M(QPHIX_ColorMatrix *m1, QPHIX_ColorMatrix *m2, QPHIX_evenodd_t parity);

void QPHIX_M_eq_V_times_Va(QPHIX_ColorMatrix *m1, QPHIX_ColorVector *v1, QPHIX_ColorVector *v2);

void QPHIX_M_meq_M(QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2);

void QPHIX_M_meq_M_times_M(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M2, QPHIX_ColorMatrix *M1);

void QPHIX_M_eq_Ma_times_Ma(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2);


/* ****************** reunit primitives ************** */
#define QPHIX_c_eq_r_plus_ir(c,a,b) {QPHIX_real(c) = (a); QPHIX_imag(c) = (b);}
#define QPHIX_c_eq_c_times_c(c,a,b)  {QPHIX_real(c) = QPHIX_real(a)*QPHIX_real(b)\
                                   - QPHIX_imag(a)*QPHIX_imag(b); \
                                    QPHIX_imag(c) = QPHIX_real(a)*QPHIX_imag(b)\
                                   + QPHIX_imag(a)*QPHIX_real(b);}

#define QPHIX_c_peq_c(c,a)   {QPHIX_real(c) += QPHIX_real(a);\
                            QPHIX_imag(c) += QPHIX_imag(a);}

#define QPHIX_c_eq_c_minus_c(c,a,b) {QPHIX_real(c) = QPHIX_real(a) - QPHIX_real(b);\
                                  QPHIX_imag(c) = QPHIX_imag(a) - QPHIX_imag(b);}

#define QPHIX_c_eq_c_plus_c(c,a,b)  {QPHIX_real(c) = QPHIX_real(a) + QPHIX_real(b);\
                                     QPHIX_imag(c) = QPHIX_imag(a) + QPHIX_imag(b);}

#define QPHIX_c_peq_r(c,a)    QPHIX_real(c) += (a);
#define QPHIX_c_eq_r(c,a)    {QPHIX_real(c) = a;\
                            QPHIX_imag(c) = 0.;}

#define QPHIX_c_peq_c_times_r(c,a,b)   {QPHIX_real(c) += (b) * QPHIX_real(a);\
                                      QPHIX_imag(c) += (b) * QPHIX_imag(a); }

#define QPHIX_c_peq_r_times_c(c,a,b)   QPHIX_c_peq_c_times_r(c,b,a)
#define QPHIX_c_eq_r_plus_ir(c,a,b) {QPHIX_real(c) = (a); QPHIX_imag(c) = (b);}

#define QPHIX_c_eq_ca(c,a)   {QPHIX_real(c) =  QPHIX_real(a);\
                              QPHIX_imag(c) = -QPHIX_imag(a);}

#define QPHIX_c_eq_r(c,a)    {QPHIX_real(c) = a;\
                              QPHIX_imag(c) = 0.;}
#define QPHIX_c_peq_c_times_c(c,a,b) {QPHIX_real(c) += (QPHIX_real(a)*QPHIX_real(b)\
                              - QPHIX_imag(a)*QPHIX_imag(b)); \
                              QPHIX_imag(c) += (QPHIX_real(a)*QPHIX_imag(b)\
                              + QPHIX_imag(a)*QPHIX_real(b));}
#define QPHIX_c_eq_c(c,a)    {QPHIX_real(c) =  QPHIX_real(a);\
                              QPHIX_imag(c) =  QPHIX_imag(a);}

#define QPHIX_c_peq_ca_times_c(c,a,b) {QPHIX_real(c) += QPHIX_real(a)*QPHIX_real(b)\
                              + QPHIX_imag(a)*QPHIX_imag(b); \
                              QPHIX_imag(c) += QPHIX_real(a)*QPHIX_imag(b)\
                              - QPHIX_imag(a)*QPHIX_real(b);}

#define QPHIX_c_eq_c_times_r(c,a,b)  {QPHIX_real(c) = (b) * QPHIX_real(a);\
                                   QPHIX_imag(c) = (b) * QPHIX_imag(a); }

#define QPHIX_c_eq_r_times_c(c,a,b)  QPHIX_c_eq_c_times_r(c,b,a)

#define QPHIX_c_eq_r_times_c_plus_c(c,a,x,b) \
                        {QPHIX_real(c) = (a)*QPHIX_real(x) + QPHIX_real(b); \
                        QPHIX_imag(c) = (a)*QPHIX_imag(x) + QPHIX_imag(b);}
#define QPHIX_r_peq_Re_c(c,a) c += QPHIX_real(a);
#define QPHIX_c_eq_r(c,a)    {QPHIX_real(c) = a; QPHIX_imag(c) = 0.;}
#define QPHIX_c_peq_c_times_ca(c,a,b) {QPHIX_real(c) += QPHIX_real(a)*QPHIX_real(b)\
                              + QPHIX_imag(a)*QPHIX_imag(b); \
                              QPHIX_imag(c) += QPHIX_imag(a)*QPHIX_real(b)\
                              - QPHIX_real(a)*QPHIX_imag(b);}
#define QPHIX_c_eq_c(c,a)    {QPHIX_real(c) =  QPHIX_real(a);\
                            QPHIX_imag(c) =  QPHIX_imag(a);}
}
#endif
