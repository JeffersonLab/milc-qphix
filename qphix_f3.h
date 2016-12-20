/*************************************/
/* Function declarations of F-3 type */
/*************************************/

#ifndef _QPHIX_F3
#define _QPHIX_F3

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QPHIX_F3_ColorVector_struct   QPHIX_F3_ColorVector;
typedef struct QPHIX_F3_FermionLinksAsqtad_struct  QPHIX_F3_FermionLinksAsqtad;
typedef struct QPHIX_F3_GaugeField_struct   QPHIX_F3_GaugeField;
typedef struct QPHIX_F3_Force_struct   QPHIX_F3_Force;

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

// Create, extract, destroy QPhiX momentum from/to raw 
QPHIX_F3_Force *QPHIX_F3_create_F_from_raw(QPHIX_F_Real *fwdrawmom, QPHIX_evenodd_t evenodd);

void QPHIX_F3_extract_F_to_raw(QPHIX_F_Real *rawdest, QPHIX_F3_Force *src, QPHIX_evenodd_t evenodd); /* from QPhiX to raw (no allocation) */

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

// create asqtad fermion links from raw
QPHIX_F3_FermionLinksAsqtad  *QPHIX_F3_asqtad_create_L_from_raw( QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_F_Real *fatback, QPHIX_F_Real *lngback, QPHIX_evenodd_t evenodd );

// copy asqtad links to raw
void  QPHIX_F3_asqtad_extract_L_to_raw( QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_F3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd );

/*
//convert asqtad links from raw
QPHIX_F3_FermionLinksAsqtad  *QPHIX_F3_asqtad_convert_L_from_raw( QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_evenodd_t evenodd );

// convert asqtad links to raw
void  QPHIX_F3_asqtad_convert_L_to_raw( QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_F3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd );

// load asqtad links from raw
void QPHIX_F3_asqtad_load_L_from_raw( QPHIX_F3_FermionLinksAsqtad *dest, QPHIX_F_Real *fat, QPHIX_F_Real *lng, QPHIX_evenodd_t evenodd );
*/

// free color vectors
void  QPHIX_F3_destroy_V( QPHIX_F3_ColorVector *V );

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
