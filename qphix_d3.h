/*************************************/
/* Function declarations of D-3 type */
/*************************************/

#ifndef _QPHIX_D3
#define _QPHIX_D3

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QPHIX_D3_ColorMatrix_struct QPHIX_D3_ColorMatrix;
typedef struct QPHIX_D3_ColorVector_struct   QPHIX_D3_ColorVector;
typedef struct QPHIX_D3_FermionLinksAsqtad_struct  QPHIX_D3_FermionLinksAsqtad;
typedef struct QPHIX_D3_FermionLinksHisq_struct  QPHIX_D3_FermionLinksHisq;
typedef struct QPHIX_D3_GaugeField_struct   QPHIX_D3_GaugeField;
typedef struct QPHIX_D3_Force_struct   QPHIX_D3_Force;

typedef struct
{
  QPHIX_D_Real real;
  QPHIX_D_Real imag;
} QPHIX_D_Complex;

typedef struct
{ 
  QPHIX_D_Complex c[3];
} SU3_D_Vector;

typedef struct
{
  QPHIX_D_Complex e[3][3];
} SU3_D_Matrix;

typedef struct
{
  QPHIX_D_Complex m01, m02, m12;
  double m00im, m11im, m22im;
  double space;
} SU3_D_AntiHermitMat;

typedef struct
{
  QPHIX_D_Complex t4[3][3][3][3];
} QPHIX_D3_ColorTensor4;

// create an empty color matrix field
QPHIX_D3_ColorMatrix *QPHIX_D3_create_M(QPHIX_evenodd_t parity);

// Destruction of color matrix field
void QPHIX_D3_destroy_M(QPHIX_D3_ColorMatrix *mat);

// create an empty color vector field
QPHIX_D3_ColorVector *QPHIX_D3_create_V(QPHIX_evenodd_t parity);

// create color vectors from raw to qphix structure
QPHIX_D3_ColorVector  *QPHIX_D3_create_V_from_raw( QPHIX_D_Real *src, QPHIX_evenodd_t evenodd );

// copy color vectors from qphix structure to raw
void QPHIX_D3_extract_V_to_raw( QPHIX_D_Real *dest, QPHIX_D3_ColorVector  *src, QPHIX_evenodd_t eveodd );

/*
// convert color vectors from raw to qphix structure
QPHIX_D3_ColorVector  *QPHIX_D3_convert_V_from_raw( QPHIX_D_Real *src, QPHIX_evenodd_t evenodd );

// convert color vectors from qphix structure to raw
QPHIX_D_Real *QPHIX_D3_convert_V_to_raw( QPHIX_D3_ColorVector *src );
*/

// free color vectors
void  QPHIX_D3_destroy_V( QPHIX_D3_ColorVector *V );

// create asqtad fermion links from raw
QPHIX_D3_FermionLinksAsqtad  *QPHIX_D3_asqtad_create_L_from_raw( QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_D_Real *fatback, QPHIX_D_Real *lngback, QPHIX_evenodd_t evenodd );

// copy asqtad links to raw
void  QPHIX_D3_asqtad_extract_L_to_raw( QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_D3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd );

// Create, extract, destroy QPhiX momentum from/to raw 
QPHIX_D3_Force *QPHIX_D3_create_F_from_raw(QPHIX_D_Real *fwdrawmom, QPHIX_evenodd_t evenodd);

void QPHIX_D3_extract_F_to_raw(QPHIX_D_Real *rawdest, QPHIX_D3_Force *src, QPHIX_evenodd_t evenodd); /* from QPhiX to raw (no allocation) */

void QPHIX_D3_destroy_F(QPHIX_D3_Force *field);

// Create, extract, destroy QPhiX gauge field from/to raw

QPHIX_D3_GaugeField *QPHIX_D3_create_G_from_raw(QPHIX_D_Real *fwdrawgauge, QPHIX_evenodd_t evenodd);

void QPHIX_D3_extract_G_to_raw(QPHIX_D_Real *rawdest, QPHIX_D3_GaugeField *src, QPHIX_evenodd_t evenodd); /* from QPhiX to raw (no allocation) */

void QPHIX_D3_destroy_G(QPHIX_D3_GaugeField *field);

// Gauge force routine

void QPHIX_D3_symanzik_1loop_gauge_force(QPHIX_info_t *info,
					 QPHIX_D3_GaugeField *gauge,
					 QPHIX_D3_Force *force,
					 QPHIX_gauge_coeffs_t *coeffs,
					 QPHIX_D_Real eps /* time step */);

/*
//convert asqtad links from raw
QPHIX_D3_FermionLinksAsqtad  *QPHIX_D3_asqtad_convert_L_from_raw( QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_evenodd_t evenodd );

// convert asqtad links to raw
void  QPHIX_D3_asqtad_convert_L_to_raw( QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_D3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd );

// load asqtad links from raw
void QPHIX_D3_asqtad_load_L_from_raw( QPHIX_D3_FermionLinksAsqtad *dest, QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_evenodd_t evenodd );
*/

// free asqtad fermion links
void  QPHIX_D3_asqtad_destroy_L( QPHIX_D3_FermionLinksAsqtad *L );

// dslash
void QPHIX_D3_asqtad_dslash( QPHIX_D3_FermionLinksAsqtad *asqtad,
			     QPHIX_D3_ColorVector *out,
			     QPHIX_D3_ColorVector *in,
			     QPHIX_evenodd_t parity );

// inverter
int QPHIX_D3_asqtad_invert( QPHIX_info_t *info,
			    QPHIX_D3_FermionLinksAsqtad *asqtad,
			    QPHIX_invert_arg_t *inv_arg,
			    QPHIX_resid_arg_t *res_arg,
			    QPHIX_D_Real mass,
			    QPHIX_D3_ColorVector *out,
			    QPHIX_D3_ColorVector *in);

// multi-mass inverter
int QPHIX_D3_asqtad_invert_multi( QPHIX_info_t *info,
				  QPHIX_D3_FermionLinksAsqtad *asqtad,
				  QPHIX_invert_arg_t *inv_arg,
				  QPHIX_resid_arg_t *res_arg[],
				  QPHIX_D_Real *mass,
				  int num_offsets,
				  QPHIX_D3_ColorVector *out[],
				  QPHIX_D3_ColorVector *in);

// fermion force
void QPHIX_D3_hisq_force_multi(QPHIX_info_t *info, QPHIX_D3_FermionLinksHisq *flh,
                            QPHIX_D3_ColorMatrix *force[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            QPHIX_D_Real *residues, QPHIX_D3_ColorVector *x[],
                            int *n_orders_naik);

void QPHIX_D3_layout_from_su3m(QPHIX_D3_ColorMatrix *dest, SU3_D_Matrix *src);
void QPHIX_D3_layout_from_4su3m(QPHIX_D3_ColorMatrix *dest[], SU3_D_Matrix *src);
void QPHIX_D3_layout_to_su3m(SU3_D_Matrix *dest, QPHIX_D3_ColorMatrix *src);
void QPHIX_D3_layout_to_4su3m(SU3_D_Matrix *dest, QPHIX_D3_ColorMatrix *src[]);

  // Creation of HISQ fermion links
QPHIX_D3_FermionLinksHisq *
QPHIX_D3_hisq_create_L_from_4su3m(void *U_links, void *V_links, void *W_unitlinks, void *Y_unitlinks,
				  QPHIX_evenodd_t parity);

// Destruction of HISQ links
void QPHIX_D3_hisq_destroy_L(QPHIX_D3_FermionLinksHisq *L);

  // Creation of a fermion force matrix (same format as color matrix)
QPHIX_D3_ColorMatrix **
QPHIX_D3_create_F_from_anti_hermitmat( void *ptr, int SZ );

  //void QPHIX_D3_layout_from_anti_hermitmat(QPHIX_D3_ColorMatrix *dest[], void *ptr, int SZ);

  //void QPHIX_D3_layout_to_anti_hermitmat(QPHIX_D3_ColorMatrix *, anti_hermitmat *);
void QPHIX_D3_layout_to_anti_hermitmat(void *ptr, QPHIX_D3_ColorMatrix *src[], int SZ);

  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QPHIX_Precision == 'D'
#include <qphix_d3_generic.h>
#endif

#ifdef __cplusplus
};
#endif

#endif // _QPHIX_D3
