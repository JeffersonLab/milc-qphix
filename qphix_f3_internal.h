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

void reconstruct_gauge_third_row( QPHIX_F_Real *g, int t, int x, int y );

#endif /* _QPHIX_F3_INTERNAL_H */
