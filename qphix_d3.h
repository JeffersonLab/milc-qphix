/*************************************/
/* Function declarations of D-3 type */
/*************************************/

#ifndef _QPHIX_D3
#define _QPHIX_D3

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QPHIX_D3_ColorVector_struct   QPHIX_D3_ColorVector;
typedef struct QPHIX_D3_FermionLinksAsqtad_struct  QPHIX_D3_FermionLinksAsqtad;
typedef struct QPHIX_D3_GaugeField_struct   QPHIX_D3_GaugeField;
typedef struct QPHIX_D3_Force_struct   QPHIX_D3_Force;

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

// create asqtad fermion links from raw
QPHIX_D3_FermionLinksAsqtad  *QPHIX_D3_asqtad_create_L_from_raw( QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_D_Real *fatback, QPHIX_D_Real *lngback, QPHIX_evenodd_t evenodd );

// copy asqtad links to raw
void  QPHIX_D3_asqtad_extract_L_to_raw( QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_D3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd );

// Create, extract, destroy QPhiX momentum from/to raw 
QPHIX_D3_Force *QPHIX_D3_create_F_from_raw(QPHIX_D_Real *fwdrawmom, QPHIX_D_Real *bckrawmom, QPHIX_evenodd_t evenodd);

void QPHIX_D3_extract_F_to_raw(QPHIX_D_Real *rawdest, QPHIX_D3_Force *src, QPHIX_evenodd_t evenodd); /* from QPhiX to raw (no allocation) */

void QPHIX_D3_destroy_F(QPHIX_D3_Force *field);

// Create, extract, destroy QPhiX gauge field from/to raw

QPHIX_D3_GaugeField *QPHIX_D3_create_G_from_raw(QPHIX_D_Real *fwdrawgauge, QPHIX_D_Real *bckrawgauge, QPHIX_evenodd_t evenodd);

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

// free color vectors
void  QPHIX_D3_destroy_V( QPHIX_D3_ColorVector *V );

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

// gauge force

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
