#ifndef _QPHIX_F3_INTERNAL_H
#define _QPHIX_F3_INTERNAL_H

/* Define the single-precision opaque objects */
/* This file is included in qphix_internal.h */

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

  /* Asqtad datatypes */

struct QPHIX_F3_FermionLinksAsqtad_struct {
  void *lngeven;
  void *lngodd;
  void *fateven;
  void *fatodd;
};


#endif /* _QPHIX_F3_INTERNAL_H */
