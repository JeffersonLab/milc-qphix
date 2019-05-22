#ifndef LAYOUT_H
#define LAYOUT_H

extern "C"
{
#include "fermion_force.h"

int QPHIX_create_V(QPHIX_ColorVector **, QPHIX_evenodd_t);
int QPHIX_create_M(QPHIX_ColorMatrix **, QPHIX_evenodd_t);
void QPHIX_destroy_M(QPHIX_ColorMatrix *);
void QPHIX_layout_from_su3m(QPHIX_ColorMatrix *, SU3_Matrix *);
void QPHIX_layout_from_4su3m(QPHIX_ColorMatrix **, SU3_Matrix *);
void QPHIX_layout_to_su3m(SU3_Matrix *, QPHIX_ColorMatrix *);
void QPHIX_layout_to_4su3m(SU3_Matrix *, QPHIX_ColorMatrix **);
void QPHIX_layout_from_anti_hermitmat(QPHIX_ColorMatrix **, void *, int SZ);
void QPHIX_layout_to_anti_hermitmat(void *, QPHIX_ColorMatrix **, int SZ);
void QPHIX_get_milc_nbr_crd(int, int, int, int, int, int *, int *, int *, int *);
void QPHIX_get_nbr_index(int, int, int, int, int *, int *, int *);
void QPHIX_V_eq_sV(QPHIX_ColorVector *, QPHIX_ColorVector *, int, int);
void QPHIX_M_eq_sM(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *, int, int);
void QPHIX_M_eq_sM_dump(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *, int, int);
void QPHIX_M_eq_sM_print(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *, int, int);
qphix_su3_matrix* QPHIX_expose_M(QPHIX_ColorMatrix* M1);
void QPHIX_reset_M(QPHIX_ColorMatrix *M2, qphix_su3_matrix* M1);

#ifdef FF_DEBUG
void QPHIX_checksum_M(QPHIX_ColorMatrix *, int);
#endif
}

#endif
