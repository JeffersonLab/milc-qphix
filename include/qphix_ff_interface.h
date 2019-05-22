#ifndef _QPHIX_FF_INTERFACE_H
#define _QPHIX_FF_INTERFACE_H

#if (QPHIX_PrecisionInt == 1)
#define QPHIX_Real float
#error QPhiX fermion force is not supported for single precision
#endif
#if (QPHIX_PrecisionInt == 2)
#define QPHIX_Real double
#endif

#define QPHIX_MAX_NAIK 20 /*verify this value */
#define QPHIX_EPS_NAIK_ZERO {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
#define QPHIX_HISQ_COEFFS_ZERO \
  ((QPHIX_hisq_coeffs_t){1, QPHIX_EPS_NAIK_ZERO, QPHIX_UNITARIZE_U3, \
   QPHIX_UNITARIZE_RATIONAL, 0,0,0,0,0, 0,0,0,0,0,0, 0,0})

// Enumerations
typedef enum
{
  QPHIX_UNITARIZE_U3 = 0,
  QPHIX_UNITARIZE_SU3 = 1,
} QPHIX_hisq_unitarize_group_t;

typedef enum
{
  QPHIX_UNITARIZE_NONE = 0,
  QPHIX_UNITARIZE_RATIONAL = 3,
  QPHIX_UNITARIZE_ANALYTIC = 5,
} QPHIX_hisq_unitarize_method_t;

// Data Structures
typedef struct 
{
  void *even;
  void *odd;
  QPHIX_evenodd_t parity;
} QPHIX_ColorMatrix;

typedef struct 
{
  int n_naiks, WeqY;
  QPHIX_ColorMatrix *U_links[4];     // gauge links
  QPHIX_ColorMatrix *V_links[4];     // Fat7 smeared
  QPHIX_ColorMatrix *Y_unitlinks[4]; // projected to U(3),
  QPHIX_ColorMatrix *W_unitlinks[4]; // projected to SU(3)
} QPHIX_FermionLinksHisq;

typedef struct
{
  int n_naiks;
  QPHIX_Real eps_naik[QPHIX_MAX_NAIK];
  //QPHIX_hisq_unitarize_group_t ugroup;
  //QPHIX_hisq_unitarize_method_t umethod;
  int ugroup;
  int umethod;
  QPHIX_Real fat7_one_link;
  QPHIX_Real fat7_three_staple;
  QPHIX_Real fat7_five_staple;
  QPHIX_Real fat7_seven_staple;
  QPHIX_Real fat7_lepage;
  QPHIX_Real asqtad_one_link;
  QPHIX_Real asqtad_three_staple;
  QPHIX_Real asqtad_five_staple;
  QPHIX_Real asqtad_seven_staple;
  QPHIX_Real asqtad_lepage;
  QPHIX_Real asqtad_naik;
  QPHIX_Real difference_one_link;
  QPHIX_Real difference_naik;
} QPHIX_hisq_coeffs_t;

// Creation & Destruction of QPHIX_ColorMatrix
int QPHIX_F3_create_M(QPHIX_ColorMatrix **, QPHIX_evenodd_t);
int QPHIX_D3_create_M(QPHIX_ColorMatrix **, QPHIX_evenodd_t);

void QPHIX_destroy_M(QPHIX_ColorMatrix *);

// Layout conversion routines
void QPHIX_F3_layout_from_su3m(QPHIX_ColorMatrix *dest, su3_matrix *src);
void QPHIX_D3_layout_from_su3m(QPHIX_ColorMatrix *dest, su3_matrix *src);

void QPHIX_F3_layout_from_anti_hermitmat(anti_hermitmat *, QPHIX_ColorMatrix *);
void QPHIX_D3_layout_from_anti_hermitmat(anti_hermitmat *, QPHIX_ColorMatrix *);

void QPHIX_F3_layout_to_anti_hermitmat(QPHIX_ColorMatrix *, anti_hermitmat *);
void QPHIX_D3_layout_to_anti_hermitmat(QPHIX_ColorMatrix *, anti_hermitmat *);

// Fermion force wrapper
void QPHIX_F3_hisq_force_multi(QPHIX_info_t *info, QPHIX_FermionLinksHisq *flh,
                            QPHIX_ColorMatrix *force[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            QPHIX_Real *residues, QPHIX_ColorVector *x[],
                            int *n_orders_naik);
void QPHIX_D3_hisq_force_multi(QPHIX_info_t *info, QPHIX_FermionLinksHisq *flh,
                            QPHIX_ColorMatrix *force[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            QPHIX_Real *residues, QPHIX_ColorVector *x[],
                            int *n_orders_naik);

#endif