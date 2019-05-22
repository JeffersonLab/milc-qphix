#ifndef FERMION_FORCE_H
#define FERMION_FORCE_H

extern "C"
{
#include "qphix.h"
#include "su3.h"

//This info comes from ks_config.h (via ks_globals.h)
//For the time being ks_globals.h is not included as it
//conflicts when this file is included in milc-qcd.
#if (QPHIX_PrecisionInt == 1)
#define VECLEN 16
#endif
#if (QPHIX_PrecisionInt == 2)
#define VECLEN 8
#endif

#define VECLENBYTES 64

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

typedef struct
{
  int x;
  int y;
  int z;
  int t;
} site;

typedef struct
{
  QPHIX_ColorVector **u;
  QPHIX_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QPHIX_eigcg_t_V;

typedef struct
{
  int dblstored, nlinks;
  QPHIX_ColorMatrix **fatlinks;
  QPHIX_ColorMatrix **longlinks;
  QPHIX_ColorMatrix **fwdlinks;
  QPHIX_ColorMatrix **bcklinks;
  QPHIX_ColorMatrix **dbllinks;
  QPHIX_eigcg_t_V eigcg;
} QPHIX_FermionLinksAsqtad_struct;

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

typedef struct
{
  char *tag;
  QPHIX_Real value;
  void *extra;
} QPHIX_opt_t;

typedef struct
{
  int inited;
  int want_deps;
  int want_aux;
  QPHIX_hisq_unitarize_method_t umethod;
  int reunit_allow_svd;
  int reunit_svd_only;
  QPHIX_Real reunit_svd_rel_error;
  QPHIX_Real reunit_svd_abs_error;
  int svd_values_info;
  int use_fat7_lepage;
} QPHIX_hisq_links_t;

typedef struct
{
  QPHIX_Real one_link;
  QPHIX_Real three_staple;
  QPHIX_Real five_staple;
  QPHIX_Real seven_staple;
  QPHIX_Real lepage;
  QPHIX_Real naik;
} QPHIX_asqtad_coeffs_t;

typedef struct
{
  int inited;
  int fnmat_src_min;
  int veclength;
  QPHIX_Real force_filter;
} QPHIX_hisq_force_t;


extern int qphix_even_sites_on_node;
extern QPHIX_hisq_links_t QPHIX_hisq_links;


void QPHIX_hisq_force_multi(QPHIX_info_t *info, QPHIX_FermionLinksHisq *flh,
                            QPHIX_ColorMatrix *force[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            QPHIX_Real *residues, QPHIX_ColorVector *x[],
                            int *n_orders_naik);

void QPHIX_get_mid(QPHIX_info_t *info, QPHIX_ColorMatrix **mid,
                   QPHIX_ColorVector *multix[], QPHIX_Real eps[],
                   QPHIX_Real scale, int nterms, int ndist);

void QPHIX_hisq_deriv_multi(QPHIX_info_t *info, QPHIX_FermionLinksHisq *flh,
                            QPHIX_ColorMatrix *deriv[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            QPHIX_Real *residues, QPHIX_ColorVector *x[],
                            int *n_orders_naik);

void QPHIX_asqtad_deriv_dbg(QPHIX_info_t *info, QPHIX_ColorMatrix *gauge[],
                        QPHIX_ColorMatrix *deriv[], QPHIX_asqtad_coeffs_t *coef,
                        QPHIX_ColorMatrix *mid_fat[],
                        QPHIX_ColorMatrix *mid_naik[]);

void QPHIX_asqtad_deriv(QPHIX_info_t *info, QPHIX_ColorMatrix *gauge[],
                        QPHIX_ColorMatrix *deriv[], QPHIX_asqtad_coeffs_t *coef,
                        QPHIX_ColorMatrix *mid_fat[],
                        QPHIX_ColorMatrix *mid_naik[]);

void QPHIX_hisq_force_multi_reunit(QPHIX_info_t *info, QPHIX_ColorMatrix *V[4],
                                   QPHIX_ColorMatrix *Force[4],
                                   QPHIX_ColorMatrix *Force_old[4]);

int QPHIX_svd2x2bidiag(QPHIX_info_t *info, QPHIX_Real *a00, QPHIX_Real *a01,
                       QPHIX_Real *a11, QPHIX_Real U2[2][2],
                       QPHIX_Real V2[2][2]);

int QPHIX_svd3x3(QPHIX_info_t *, qphix_su3_matrix *, QPHIX_Real *,
                 qphix_su3_matrix *, qphix_su3_matrix *);

void u3_un_der_analytic(QPHIX_info_t *, qphix_su3_matrix *,
                        QPHIX_ColorTensor4 *, QPHIX_ColorTensor4 *);

QPHIX_Complex su3_mat_det(qphix_su3_matrix *);

}

#endif
