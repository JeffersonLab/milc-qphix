/*************************************/
/* Function declarations of F-3 type */
/*************************************/

#ifndef _QPHIX_F3
#define _QPHIX_F3

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QPHIX_F3_ColorMatrix_struct QPHIX_F3_ColorMatrix;
typedef struct QPHIX_F3_ColorVector_struct   QPHIX_F3_ColorVector;
typedef struct QPHIX_F3_FermionLinksAsqtad_struct  QPHIX_F3_FermionLinksAsqtad;
typedef struct QPHIX_F3_FermionLinksHisq_struct  QPHIX_F3_FermionLinksHisq;
typedef struct QPHIX_F3_GaugeField_struct   QPHIX_F3_GaugeField;
typedef struct QPHIX_F3_Force_struct   QPHIX_F3_Force;

typedef struct
{
  QPHIX_F_Real real;
  QPHIX_F_Real imag;
} QPHIX_F_Complex;

typedef struct
{ 
  QPHIX_F_Complex c[3];
} SU3_F_Vector;

typedef struct
{
  QPHIX_F_Complex e[3][3];
} SU3_F_Matrix;

typedef struct
{
  QPHIX_F_Complex m01, m02, m12;
  float m00im, m11im, m22im;
  float space;
} SU3_F_AntiHermitMat;

typedef struct
{
  QPHIX_F_Complex t4[3][3][3][3];
} QPHIX_F3_ColorTensor4;

#include "veclength.h"
typedef QPHIX_F_Real QSU3V_float[3][2][VECLEN];
typedef QPHIX_F_Real QSU3M_float[3][3][2][VECLEN];
typedef QPHIX_F_Complex qphix_float_su3_matrix[3][3];

// create an empty color matrix field
QPHIX_F3_ColorMatrix *QPHIX_F3_create_M(QPHIX_evenodd_t parity);

// Destruction of color matrix field
void QPHIX_F3_destroy_M(QPHIX_F3_ColorMatrix *mat);

// create color vectors from raw to qphix structure
QPHIX_F3_ColorVector  *QPHIX_F3_create_V_from_raw( QPHIX_F_Real *src, QPHIX_evenodd_t evenodd );

// copy color vectors from qphix structure to raw
void QPHIX_F3_extract_V_to_raw( QPHIX_F_Real *dest, QPHIX_F3_ColorVector  *src, QPHIX_evenodd_t evenodd );

/*
// convert color vectors from raw to qphix structure
QPHIX_F3_ColorVector  *QPHIX_F3_convert_V_from_raw( QPHIX_F_Real *src, QPHIX_evenodd_t evenodd );

// convert color vectors from qphix structure to raw
QPHIX_F_Real *QPHIX_F3_convert_V_to_raw( QPHIX_F3_ColorVector *src );
*/

// free color vectors
void  QPHIX_F3_destroy_V( QPHIX_F3_ColorVector *V );

// create asqtad fermion links from raw
QPHIX_F3_FermionLinksAsqtad  *QPHIX_F3_asqtad_create_L_from_raw( QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_F_Real *fatback, QPHIX_F_Real *lngback, QPHIX_evenodd_t evenodd );

// copy asqtad links to raw
void  QPHIX_F3_asqtad_extract_L_to_raw( QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_F3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd );

// Create, extract, destroy QPhiX momentum from/to raw 
QPHIX_F3_Force *QPHIX_F3_create_F_from_raw(QPHIX_F_Real *fwdrawmom, QPHIX_evenodd_t evenodd);

void QPHIX_F3_extract_F_to_raw(QPHIX_F_Real *rawdest, QPHIX_F3_Force *src, QPHIX_evenodd_t evenodd); /* from QPhiX to raw (no allocation) */

// Free fermion force 
void QPHIX_F3_destroy_F(QPHIX_F3_Force *field);

// Create, extract, destroy QPhiX gauge field from/to raw

QPHIX_F3_GaugeField *QPHIX_F3_create_G_from_raw(QPHIX_F_Real *fwdrawgauge, QPHIX_evenodd_t evenodd);

void QPHIX_F3_extract_G_to_raw(QPHIX_F_Real *rawdest, QPHIX_F3_GaugeField *src, QPHIX_evenodd_t evenodd); /* from QPhiX to raw (no allocation) */

void QPHIX_F3_destroy_G(QPHIX_F3_GaugeField *field);

// Gauge force routine

void QPHIX_F3_symanzik_1loop_gauge_force(QPHIX_info_t *info,
					 QPHIX_F3_GaugeField *gauge,
					 QPHIX_F3_Force *force,
					 QPHIX_gauge_coeffs_t *coeffs,
					 QPHIX_F_Real eps /* time step */);

/*
//convert asqtad links from raw
QPHIX_F3_FermionLinksAsqtad  *QPHIX_F3_asqtad_convert_L_from_raw( QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_evenodd_t evenodd );

// convert asqtad links to raw
void  QPHIX_F3_asqtad_convert_L_to_raw( QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_F3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd );

// load asqtad links from raw
void QPHIX_F3_asqtad_load_L_from_raw( QPHIX_F3_FermionLinksAsqtad *dest, QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_evenodd_t evenodd );
*/

// free asqtad fermion links
void  QPHIX_F3_asqtad_destroy_L( QPHIX_F3_FermionLinksAsqtad *L );

// dslash
void QPHIX_F3_asqtad_dslash( QPHIX_F3_FermionLinksAsqtad *asqtad,
			     QPHIX_F3_ColorVector *out,
			     QPHIX_F3_ColorVector *in,
			     QPHIX_evenodd_t parity );

// inverter
int QPHIX_F3_asqtad_invert( QPHIX_info_t *info,
			    QPHIX_F3_FermionLinksAsqtad *asqtad,
			    QPHIX_invert_arg_t *inv_arg,
			    QPHIX_resid_arg_t *res_arg,
			    QPHIX_F_Real mass,
			    QPHIX_F3_ColorVector *out,
			    QPHIX_F3_ColorVector *in);

// multi-mass inverter
int QPHIX_F3_asqtad_invert_multi( QPHIX_info_t *info,
				  QPHIX_F3_FermionLinksAsqtad *asqtad,
				  QPHIX_invert_arg_t *inv_arg,
				  QPHIX_resid_arg_t *res_arg[],
				  QPHIX_F_Real *mass,
				  int num_offsets,
				  QPHIX_F3_ColorVector *out[],
				  QPHIX_F3_ColorVector *in);
// fermion force
void QPHIX_F3_hisq_force_multi(QPHIX_info_t *info, QPHIX_F3_FermionLinksHisq *flh,
                            QPHIX_F3_ColorMatrix *force[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            QPHIX_F_Real *residues, QPHIX_F3_ColorVector *x[],
                            int *n_orders_naik);


#include "su3.h"
void QPHIX_F3_layout_from_su3m(QPHIX_F3_ColorMatrix *dest, SU3_F_Matrix *src);
void QPHIX_F3_layout_from_4su3m(QPHIX_F3_ColorMatrix *dest[], SU3_F_Matrix *src);
void QPHIX_F3_layout_to_su3m(SU3_F_Matrix *dest, QPHIX_F3_ColorMatrix *src);
void QPHIX_F3_layout_to_4su3m(SU3_F_Matrix *dest, QPHIX_F3_ColorMatrix *src[]);


  // Creation of HISQ fermion links
QPHIX_F3_FermionLinksHisq *
QPHIX_F3_hisq_create_L_from_4su3m(void *U_links, void *V_links, void *W_unitlinks, void *Y_unitlinks,
				  QPHIX_evenodd_t parity);

// Destruction of HISQ links
void QPHIX_F3_hisq_destroy_L(QPHIX_F3_FermionLinksHisq *L);

  // Creation of a fermion force matrix (same format as color matrix)
QPHIX_F3_ColorMatrix **
QPHIX_F3_create_F_from_anti_hermitmat( void *ptr, int SZ );

  //void QPHIX_F3_layout_from_anti_hermitmat(anti_hermitmat *, QPHIX_F3_ColorMatrix *);
  //void QPHIX_F3_layout_from_anti_hermitmat(QPHIX_F3_ColorMatrix *dest[], void *ptr, int SZ);

  //void QPHIX_F3_layout_to_anti_hermitmat(QPHIX_F3_ColorMatrix *, anti_hermitmat *);
void QPHIX_F3_layout_to_anti_hermitmat(void *ptr, QPHIX_F3_ColorMatrix *src[], int SZ);
void QPHIX_F3_reset_M(QPHIX_F3_ColorMatrix *M2, qphix_float_su3_matrix* M1);
qphix_float_su3_matrix* QPHIX_F3_expose_M(QPHIX_F3_ColorMatrix* M1);

  // Internal functions for the fermion force
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



  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QPHIX_Precision == 'F'
#include <qphix_f3_generic.h>
#endif

#ifdef __cplusplus
};
#endif

#endif // _QPHIX_F3
