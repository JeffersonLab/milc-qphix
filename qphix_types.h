#ifndef _QPHIX_TYPES
#define _QPHIX_TYPES

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct QPHIX_F3_ColorVector {
    void *even;
    void *odd;
    int parity;
  } QPHIX_F3_ColorVector;

  typedef struct QPHIX_F3_FermionLinksAsqtad {
    void *lngeven;
    void *lngodd;
    void *fateven;
    void *fatodd;
  } QPHIX_F3_FermionLinksAsqtad;

  typedef float QPHIX_F_Real;
  typedef double QPHIX_D_Real;
  typedef int QPHIX_evenodd_t;

#ifdef __cplusplus
};
#endif

#endif
