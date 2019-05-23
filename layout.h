#ifndef LAYOUT_H
#define LAYOUT_H

extern "C"
{
#include "fermion_force.h"

void QPHIX_get_milc_nbr_crd(int, int, int, int, int, int *, int *, int *, int *);
void QPHIX_get_nbr_index(int, int, int, int, int *, int *, int *);
void QPHIX_V_eq_sV(QPHIX_ColorVector *, QPHIX_ColorVector *, int, int);
void QPHIX_M_eq_sM(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *, int, int);
void QPHIX_M_eq_sM_dump(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *, int, int);
void QPHIX_M_eq_sM_print(QPHIX_ColorMatrix *, QPHIX_ColorMatrix *, int, int);

#ifdef FF_DEBUG
void QPHIX_checksum_M(QPHIX_ColorMatrix *, int);
#endif
}

#endif
