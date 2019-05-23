#ifndef _QPHIX_F3_INTERNAL_H
#define _QPHIX_F3_INTERNAL_H

/* Define the single-precision opaque objects */
/* This file is included in qphix_internal.h */

struct QPHIX_F3_ColorMatrix_struct {
  void *even;
  void *odd;
  int parity;
};

struct QPHIX_F3_ColorVector_struct {
  void *even;
  void *odd;
  int parity;
};

struct QPHIX_F3_GaugeField_struct {
  void *geo;
  int parity;
  int compress12;
};

struct QPHIX_F3_Force_struct {
  void *heo;
  void *htmp;
  int parity;
};

  /* Asqtad datatypes used by the inverters */

struct QPHIX_F3_FermionLinksAsqtad_struct {
  void *lngeven;
  void *lngodd;
  void *fateven;
  void *fatodd;
};

  /* HISQ datatypes */

struct  QPHIX_F3_FermionLinksHisq_struct{
  int n_naiks, WeqY;
  QPHIX_F3_ColorMatrix *U_links[4];     // gauge links
  QPHIX_F3_ColorMatrix *V_links[4];     // Fat7 smeared
  QPHIX_F3_ColorMatrix *Y_unitlinks[4]; // projected to U(3),
  QPHIX_F3_ColorMatrix *W_unitlinks[4]; // projected to SU(3)
};

/* Internal functions for the fermion force */

typedef QPHIX_F_Real QSU3V_float[3][2][VECLEN];
typedef QPHIX_F_Real QSU3M_float[3][3][2][VECLEN];
typedef QPHIX_F_Complex qphix_float_su3_matrix[3][3];

#if (QPHIX_PrecisionInt == 1)
#define QSU3M QSU3M_float
#define QSU3V QSU3V_float
#define qphix_su3_matrix qphix_float_su3_matrix
#endif

void QPHIX_F3_get_mid(QPHIX_info_t *info, QPHIX_F3_ColorMatrix **mid,
                   QPHIX_F3_ColorVector *multix[], float eps[],
                   float scale, int nterms, int ndist);

void QPHIX_F3_hisq_deriv_multi(QPHIX_info_t *info, QPHIX_F3_FermionLinksHisq *flh,
                            QPHIX_F3_ColorMatrix *deriv[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            float *residues, QPHIX_F3_ColorVector *x[],
                            int *n_orders_naik);

void QPHIX_F3_asqtad_deriv_dbg(QPHIX_info_t *info, QPHIX_F3_ColorMatrix *gauge[],
                        QPHIX_F3_ColorMatrix *deriv[], QPHIX_asqtad_coeffs_t *coef,
                        QPHIX_F3_ColorMatrix *mid_fat[],
                        QPHIX_F3_ColorMatrix *mid_naik[]);

void QPHIX_F3_asqtad_deriv(QPHIX_info_t *info, QPHIX_F3_ColorMatrix *gauge[],
                        QPHIX_F3_ColorMatrix *deriv[], QPHIX_asqtad_coeffs_t *coef,
                        QPHIX_F3_ColorMatrix *mid_fat[],
                        QPHIX_F3_ColorMatrix *mid_naik[]);

void QPHIX_F3_hisq_force_multi_reunit(QPHIX_info_t *info, QPHIX_F3_ColorMatrix *V[4],
                                   QPHIX_F3_ColorMatrix *Force[4],
                                   QPHIX_F3_ColorMatrix *Force_old[4]);

int QPHIX_F3_svd2x2bidiag(QPHIX_info_t *info, float *a00, float *a01,
                       float *a11, float U2[2][2],
                       float V2[2][2]);

int QPHIX_F3_svd3x3(QPHIX_info_t *, qphix_float_su3_matrix *, float *,
                 qphix_float_su3_matrix *, qphix_float_su3_matrix *);

void u3_un_der_analytic_float(QPHIX_info_t *, qphix_float_su3_matrix *,
			       QPHIX_F3_ColorTensor4 *, QPHIX_F3_ColorTensor4 *);

QPHIX_F_Complex su3_mat_det_float(qphix_float_su3_matrix *);

void reconstruct_gauge_third_row( QPHIX_F_Real *g, int t, int x, int y );

#endif /* _QPHIX_F3_INTERNAL_H */
