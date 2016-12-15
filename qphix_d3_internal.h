#ifndef _QPHIX_D3_INTERNAL_H
#define _QPHIX_D3_INTERNAL_H

/* Define the single-precision opaque objects */
/* This file is included in qphix_internal.h */

struct QPHIX_D3_ColorVector_struct {
  void *even;
  void *odd;
  int parity;
};

  /* Asqtad datatypes */

struct QPHIX_D3_FermionLinksAsqtad_struct {
  void *lngeven;
  void *lngodd;
  void *fateven;
  void *fatodd;
};

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

#endif /* _QPHIX_D3_INTERNAL_H */
