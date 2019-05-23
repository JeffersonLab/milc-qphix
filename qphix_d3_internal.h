#ifndef _QPHIX_D3_INTERNAL_H
#define _QPHIX_D3_INTERNAL_H

/* Define the single-precision opaque objects */
/* This file is included in qphix_internal.h */

struct QPHIX_D3_ColorMatrix_struct {
  void *even;
  void *odd;
  int parity;
};

struct QPHIX_D3_ColorVector_struct {
  void *even;
  void *odd;
  int parity;
};

  /* Asqtad datatypes used by the inverters */

struct QPHIX_D3_FermionLinksAsqtad_struct {
  void *lngeven;
  void *lngodd;
  void *fateven;
  void *fatodd;
};

  /* HISQ datatypes */

struct  QPHIX_D3_FermionLinksHisq_struct{
  int n_naiks, WeqY;
  QPHIX_D3_ColorMatrix *U_links[4];     // gauge links
  QPHIX_D3_ColorMatrix *V_links[4];     // Fat7 smeared
  QPHIX_D3_ColorMatrix *Y_unitlinks[4]; // projected to U(3),
  QPHIX_D3_ColorMatrix *W_unitlinks[4]; // projected to SU(3)
};

  /* Gauge field */

struct QPHIX_D3_GaugeField_struct {
  void *geo;
  int parity;
  int compress12;
};

struct QPHIX_D3_Force_struct {
  void *heo;
  void *htmp;
  int parity;
};

/* Internal functions for the fermion force */

typedef QPHIX_D_Real QSU3V_double[3][2][VECLEN];
typedef QPHIX_D_Real QSU3M_double[3][3][2][VECLEN];
typedef QPHIX_D_Complex qphix_double_su3_matrix[3][3];

#if (QPHIX_PrecisionInt == 2)
#define QSU3M QSU3M_double
#define QSU3V QSU3V_double
#define qphix_su3_matrix qphix_double_su3_matrix
#endif

void QPHIX_D3_get_mid(QPHIX_info_t *info, QPHIX_D3_ColorMatrix **mid,
                   QPHIX_D3_ColorVector *multix[], double eps[],
                   double scale, int nterms, int ndist);

void QPHIX_D3_hisq_deriv_multi(QPHIX_info_t *info, QPHIX_D3_FermionLinksHisq *flh,
                            QPHIX_D3_ColorMatrix *deriv[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            double *residues, QPHIX_D3_ColorVector *x[],
                            int *n_orders_naik);

void QPHIX_D3_asqtad_deriv_dbg(QPHIX_info_t *info, QPHIX_D3_ColorMatrix *gauge[],
                        QPHIX_D3_ColorMatrix *deriv[], QPHIX_asqtad_coeffs_t *coef,
                        QPHIX_D3_ColorMatrix *mid_fat[],
                        QPHIX_D3_ColorMatrix *mid_naik[]);

void QPHIX_D3_asqtad_deriv(QPHIX_info_t *info, QPHIX_D3_ColorMatrix *gauge[],
                        QPHIX_D3_ColorMatrix *deriv[], QPHIX_asqtad_coeffs_t *coef,
                        QPHIX_D3_ColorMatrix *mid_fat[],
                        QPHIX_D3_ColorMatrix *mid_naik[]);

void QPHIX_D3_hisq_force_multi_reunit(QPHIX_info_t *info, QPHIX_D3_ColorMatrix *V[4],
                                   QPHIX_D3_ColorMatrix *Force[4],
                                   QPHIX_D3_ColorMatrix *Force_old[4]);

int QPHIX_D3_svd2x2bidiag(QPHIX_info_t *info, double *a00, double *a01,
                       double *a11, double U2[2][2],
                       double V2[2][2]);

int QPHIX_D3_svd3x3(QPHIX_info_t *, qphix_double_su3_matrix *, double *,
                 qphix_double_su3_matrix *, qphix_double_su3_matrix *);

void u3_un_der_analytic_double(QPHIX_info_t *, qphix_double_su3_matrix *,
			       QPHIX_D3_ColorTensor4 *, QPHIX_D3_ColorTensor4 *);

QPHIX_D_Complex su3_mat_det_double(qphix_double_su3_matrix *);

void reconstruct_gauge_third_row( QPHIX_D_Real *g, int t, int x, int y );

#endif /* _QPHIX_D3_INTERNAL_H */
