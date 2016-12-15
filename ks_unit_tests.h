/****************************** Unit Test functions **********************/
#ifndef _KS_UNIT_TESTS_H_
#define _KS_UNIT_TESTS_H_

#ifdef __cplusplus
extern "C" {
#endif

double
t_magsq_su3vec (void *ks_spinor);

void 
t_scalar_mult_add_su3_vector (const void *sp1, const void *sp2 
                              , double s, void *sp3);

void 
t_scalar_mult_su3_vector (const void *sp1, double s, void *sp2);

double
t_su3_rdot (void *sp1, void *sp2);

void 
t_add_su3_vector (const void *sp1, const void *sp2, void *sp3);

void 
t_sub_su3_vector (const void *sp1, const void *sp2, void *sp3);


double 
t_relative_residue (void *p, void *q, int parity);

#ifdef __cplusplus
}
#endif

#endif
