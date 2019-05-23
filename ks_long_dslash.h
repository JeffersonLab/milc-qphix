#ifndef _KS_LONG_DSLASH_H_
#define _KS_LONG_DSLASH_H_
/*
#define EVEN 0x02
#define ODD 0x01
#define EVENANDODD 0x03

#define OPP_PAR(p) ((p)==EVENANDODD ? (p) : (p)^11)
*/
#ifdef __cplusplus
extern "C" {
#endif

#include "ks_unit_tests.h"

extern const int MILC2MBENCH[];
extern const int MBENCH2MILC[];

typedef struct {
	unsigned char x, y, z, t, d;
} ind_t;

void 
dump_ks_spinor (void * spinor);

/*
void 
setup_mbench( int *lattice, int *geom, int *my_coord, int *my_neigh_ranks, int my_id
             , int NCores_, int ThreadsPerCore, int minCt_);
*/

void 
*allocKS_F ();

void 
*allocKS2_F ();

void freeKS (void *p);

void 
*allocGauge_F ();

void 
*allocGauge18_F ();

void 
freeGauge (void *p);

void 
freeGauge18 (void *p);

void 
setKS_F (void *in, void *su3_v, int x, int y, int z, int t);

void 
getKS_F (void *out, void *su3_v, int x, int y, int z, int t);

void 
setGauge_F (void *in, void *su3_m, int dir, int x, int y, int z, int t);

void 
getGauge_F (void *out, void *su3_m, int dir, int x, int y, int z, int t);

void 
setGauge18_F (void *in, void *su3_m, int dir, int x, int y, int z, int t);

void 
getGauge18_F (void *out, void *su3_m, int dir, int x, int y, int z, int t);

void 
setFullGauge_F (void *in, void *su3_m, int len);

void 
setFullGauge18_F (void *in, void *su3_m, int len);

void 
qphix_ks_dslash_F (void *s_in, void *gll_, void *gfl_, void *s_out_, int cb);

void
*allocKS_F ();

void
*allocKS2_F ();

void freeKS (void *p);

void
*allocGauge_F ();

void
*allocGauge18_F ();

void
freeGauge (void *p);

void
freeGauge18 (void *p);

void
setKS_F (void *in, void *su3_v, int x, int y, int z, int t);

void
getKS_F (void *out, void *su3_v, int x, int y, int z, int t);

void
setGauge_F (void *in, void *su3_m, int dir, int x, int y, int z, int t);

void
getGauge_F (void *out, void *su3_m, int dir, int x, int y, int z, int t);

void
setGauge18_F (void *in, void *su3_m, int dir, int x, int y, int z, int t);

void
getGauge18_F (void *out, void *su3_m, int dir, int x, int y, int z, int t);

void
setFullGauge_F (void *in, void *su3_m, int len);

void
setFullGauge18_F (void *in, void *su3_m, int len);

void
qphix_ks_dslash_F (void *s_in, void *gll_, void *gfl_, void *s_out_, int cb);

void
*allocKS_D ();

void
*allocKS2_D ();

void
*allocGauge_D ();

void
*allocGauge18_D ();

void
setKS_D (void *in, void *su3_v, int x, int y, int z, int t);

void
getKS_D (void *out, void *su3_v, int x, int y, int z, int t);

void
setGauge_D (void *in, void *su3_m, int dir, int x, int y, int z, int t);

void
getGauge_D (void *out, void *su3_m, int dir, int x, int y, int z, int t);

void
setGauge18_D (void *in, void *su3_m, int dir, int x, int y, int z, int t);


void
getGauge18_D (void *out, void *su3_m, int dir, int x, int y, int z, int t);

void
setFullGauge_D (void *in, void *su3_m, int len);

void
setFullGauge18_D (void *in, void *su3_m, int len);

void
qphix_ks_dslash_D (void *s_in, void *gll_, void *gfl_, void *s_out_, int cb);

#ifdef __cplusplus
}
#endif

#if QPHIX_PrecisionInt==1
#define allocKS allocKS_F
#define allocKS2 allocKS2_F
#define allocGauge allocGauge_F
#define allocGauge18 allocGauge18_F
#define setKS setKS_F
#define getKS getKS_F
#define setFullGauge18 setFullGauge18_F
#define setFullGauge setFullGauge_F
#define setGauge setGauge_F
#define getGauge getGauge_F
#define setGauge18 setGauge18_F
#define getGauge18 getGauge18_F
#define qphix_ks_dslash qphix_ks_dslash_F
#elif QPHIX_PrecisionInt==2
#define allocKS allocKS_D
#define allocKS2 allocKS2_D
#define allocGauge allocGauge_D
#define allocGauge18 allocGauge18_D
#define setKS setKS_D
#define getKS getKS_D
#define setFullGauge18 setFullGauge18_D
#define setFullGauge setFullGauge_D
#define setGauge setGauge_D
#define getGauge getGauge_D
#define setGauge18 setGauge18_D
#define getGauge18 getGauge18_D
#define qphix_ks_dslash qphix_ks_dslash_D
#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif

#endif // _KS_LONG_DSLASH_H_
