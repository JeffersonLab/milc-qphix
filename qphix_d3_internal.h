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

void reconstruct_gauge_third_row( QPHIX_D_Real *g, int t, int x, int y );

#endif /* _QPHIX_D3_INTERNAL_H */
