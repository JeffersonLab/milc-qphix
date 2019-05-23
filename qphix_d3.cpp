/*************************************/
/* Definitions of F-3 type functions */
/*************************************/
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <immintrin.h>
#include <iostream>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ks_config.h"
#include "ks_globals.h"        /* All global values and macros */
#include "ks_long_dslash.h"
#include "gauge_force_imp.h"
#include "ks_boundary.h"
#include "fermion_force.h"

#include "qphix_internal.h"
#define TIME_CG 1
#define TIME_GF 1
#include <sys/time.h>
#if TIME_CG
#include "hrtimer.h"   /* High resolution timer */
#include "ks_cg_perf_profiler.h"
#endif
#if TIME_GF
#include "hrtimer.h"
#endif

/*
#if QPHIX_PrecisionInt == 1
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
*/
#if QPHIX_PrecisionInt == 2
#define allocKS allocKS_D
#define allocKS2 allocKS2_D
#define allocGauge allocGauge_D
#define allocGauge18 allocGauge18_D
#define allocHermit allocHermit_D
#define allocHermitHelperYZT allocHermitHelperYZT_D
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

const int nGX = (VECLEN < 2 ? 1 : 2);
const int nGY = (VECLEN < 4 ? 1 : 2);
const int nGZ = (VECLEN < 8 ? 1 : 2);
const int nGT = (VECLEN < 16 ? 1 : 2);
extern int qphix_even_sites_on_node;
extern int qphix_fused_sites_on_node;
//extern int sites_on_node;

static hrtimer_t timer = 0;

/* Creation of empty QPHIX_ColorMatrix field */
QPHIX_D3_ColorMatrix *
QPHIX_D3_create_M(QPHIX_evenodd_t parity)
{
  int p;
  int pmap[2];
  
  pmap[0] = parity & QPHIX_EVEN;
  pmap[1] = parity & QPHIX_ODD;
  
  QPHIX_ColorMatrix *m1 = (QPHIX_ColorMatrix *) malloc(sizeof(QPHIX_ColorMatrix));
  MYASSERT(m1 != NULL)
  m1->even = NULL;
  m1->odd  = NULL;
  m1->parity = parity;
  
//  for(p=0; p<2; p++)
  for(p=0; p<1; p++)
  {
    if (pmap[0])
    {
      m1->even = (QSU3M_double *) _mm_malloc(qphix_even_sites_on_node * sizeof(QSU3M_double), VECLENBYTES);
      MYASSERT(m1->even != NULL);
    }
    if (pmap[1])
    {
      m1->odd  = (QSU3M_double *) _mm_malloc(qphix_even_sites_on_node * sizeof(QSU3M_double), VECLENBYTES);
      MYASSERT(m1->odd != NULL);
    }
  }
  return m1;
}

/* Destruction of QPHIX_ColorMatrix */
void QPHIX_D3_destroy_M(QPHIX_D3_ColorMatrix *mat)
{
  if (mat == NULL) return;
  if (mat->even != NULL) _mm_free(mat->even);
  if (mat->odd  != NULL) _mm_free(mat->odd);
  free(mat);
}

/* Creation of empty QPHIX_ColorVector field */
QPHIX_D3_ColorVector *
QPHIX_D3_create_V(QPHIX_evenodd_t parity)
{
    int p;
    int pmap[2];

    pmap[0] = parity & QPHIX_EVEN;
    pmap[1] = parity & QPHIX_ODD;
    
    QPHIX_D3_ColorVector *v1 = (QPHIX_ColorVector *) malloc(sizeof(QPHIX_ColorVector));
    MYASSERT(v1 != NULL);

    v1->even = NULL;
    v1->odd  = NULL;
    v1->parity = parity;
    
    for(p=0; p<2; p++)
    {
      if (pmap[0])
      {
        v1->even = (QSU3V *) _mm_malloc(qphix_even_sites_on_node * sizeof(QSU3V), VECLENBYTES);
        MYASSERT(v1->even != NULL);
      }
      if (pmap[1])
      {
        v1->odd  = (QSU3V *) _mm_malloc(qphix_even_sites_on_node * sizeof(QSU3V), VECLENBYTES);
        MYASSERT(v1->odd != NULL)
      }
    }
    return v1;
  }

QPHIX_D3_ColorVector *
QPHIX_D3_create_V_from_raw( QPHIX_D_Real *src, QPHIX_evenodd_t evenodd )
{
    QPHIX_D3_ColorVector *V = (QPHIX_D3_ColorVector *)malloc(sizeof(QPHIX_D3_ColorVector));

    MYASSERT(V!=0x00);
    V->parity = evenodd;
    V->even = 0x00;
    V->odd = 0x00;
    int peven = evenodd & QPHIX_EVEN;
    int podd = evenodd & QPHIX_ODD;

#if TIME_CG
    timer = -hrtimer_get();
#endif
    if( peven ) 
    {
	V->even = allocKS();
	MYASSERT(V->even!=0x00);
    }
    if( podd )
    {
	V->odd = allocKS();
	MYASSERT(V->odd!=0x00);
    }
    KS *kse = (KS *)(V->even);
    KS *kso = (KS *)(V->odd);
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads)     
#endif
    for(int i=0; i<qphix_even_sites_on_node; ++i) 
    {
	    int y = i / Nxh;
	    int x = i - y * Nxh;
	    int z = y / Ny;
	    y -= z * Ny;
	    int t = z / Nz;
	    z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
	    x1 = x / Vxh; x2 = x % Vxh;
	    y1 = y / Vy; y2 = y % Vy;
	    z1 = z / Vz; z2 = z % Vz;
	    t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

	    for(int c = 0; c < 3; c++) {
		if(peven) {
		    kse[ind][c][0][v] = src[6*i+c*2];
		    kse[ind][c][1][v] = src[6*i+c*2+1];
		}
		if(podd) {
		    kso[ind][c][0][v] = src[6*(i+qphix_even_sites_on_node)+c*2];
		    kso[ind][c][1][v] = src[6*(i+qphix_even_sites_on_node)+c*2+1];
		}
	    }
    }
/*
    if( evenodd & QPHIX_ODD ) {
	V->odd = initBuf(_mm_malloc((Pxyz*Vt+1)*sizeof(KS), 64), (Pxyz*Vt+1)*sizeof(KS));
	KS *ks = (KS *)(V->odd);
#ifdef _OPENMP
#pragma omp for default(shared) num_threads(nThreads)     
#endif
	for( int i=0; i<n_odd_sites_on_node; ++i ) {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = (x-Lsxh) / Vxh; x2 = x % Vxh;
            y1 = (y-Lsy) / Vy; y2 = y % Vy;
            z1 = (z-Lsz) / Vz; z2 = z % Vz;
            t1 = (t-Lst) / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

            for(int c = 0; c < 3; c++) {
		ks[ind][c][0][v] = src[6*i+2*c+qphix_even_sites_on_node];
		ks[ind][c][1][v] = src[6*i+2*c+1+qphix_even_sites_on_node];
	    }
	}
    }
*/
#if TIME_CG
    timer +=  hrtimer_get();
    //if(myRank==0) printf("QPHIX_D3_create_V_from_raw: time = %.6e sec\n ", (timer)/1.e+9);
    fflush(stdout);
#endif
    return V;
}

/* create SU3 gauge fields from raw to qphix structure */
QPHIX_D3_GaugeField  
*QPHIX_D3_create_G_from_raw( QPHIX_D_Real *rawsrc, QPHIX_evenodd_t evenodd )
{
    MYASSERT(evenodd==QPHIX_EVENODD||evenodd==QPHIX_EVEN||evenodd==QPHIX_ODD);
#if TIME_GF
    timer = -hrtimer_get();
#endif
    QPHIX_D3_GaugeField *gf = (QPHIX_D3_GaugeField *)malloc(sizeof(QPHIX_D3_GaugeField));
    MYASSERT(gf!=0x00);
    gf->geo = allocGauge();
    MYASSERT(gf->geo!=0x00);
#if COMPRESSED_12
    gf->compress12 = 1;
#else
    gf->compress12 = 0;
#endif
    gf->parity = evenodd;
    /* work on EVEN sites if parity is QPHIX_EVENODD */
    if(evenodd == QPHIX_EVENODD) gf->parity = QPHIX_EVEN;

    Gauge *dest = (Gauge *)gf->geo;
    fptype *rawsrcp;
    //if(evenodd & QPHIX_ODD) rawsrcp = rawsrc + 72*qphix_even_sites_on_node;
    rawsrcp = rawsrc + 72*qphix_even_sites_on_node;
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    if(gf->parity == QPHIX_EVEN)
    {
        /*if(evenodd & QPHIX_ODD)*/ 
	for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
        {
            pack_gauge_face_dir(tid, rawsrcp, 2*dir+1, 1);
        }

    	for(int i=tid; i<qphix_even_sites_on_node; i+=nThreads)
    	{
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;
	    int xodd = (y + z + t) & 1;

	    for(int d=0; d<4; ++d) 
	    {
                int backi = i;
                if(d==0) {
                    backi -= 1-xodd;
                    if(xodd==0 && x==0) backi += Nxh;
                }
                else if(d==1) {
                    backi -= Nxh;
                    if(y==0) backi += Nxh*Ny;
                }
                else if(d==2) {
                    backi -= Nxh*Ny;
                    if(z==0) backi += Nxh*Ny*Nz;
                }
                else {
                    backi -= Nxh*Ny*Nz;
                    if(t==0) backi += Nxh*Ny*Nz*Nt;
                }
#if COMPRESSED_12
	    	for(int c1=0; c1<2; ++c1) 
#else
	    	for(int c1=0; c1<3; ++c1)
#endif
	    	for(int c2=0; c2<3; ++c2) for(int ir=0; ir<2; ++ir) 
	    	{
		    dest[ind][2*d+1][c1][c2][ir][v] = rawsrc[72*i+18*d+c1*6+c2*2+ir];
		    //if(evenodd & QPHIX_ODD)
			dest[ind][2*d][c1][c2][ir][v] = rawsrcp[backi*72+18*d+c1*6+c2*2+ir];
		}
	    }
	}
	/*if(evenodd & QPHIX_ODD)*/
 	for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
        {
            unpack_gauge_face_dir(tid, dest, 2*dir, 0);
        }
    }
    else
    {
	for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
	{
	    pack_gauge_face_dir(tid, rawsrc, 2*dir+1, 0);
	}
    	for(int i=tid; i<qphix_even_sites_on_node; i+=nThreads)
    	{
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;
	    int xodd = (y + z + t) & 1;

            for(int d=0; d<4; ++d)
	    {
                int backi = i;
                if(d==0) {
                    backi -= xodd;
                    if(xodd==1 && x==0) backi += Nxh;
                }
                else if(d==1) {
                    backi -= Nxh;
                    if(y==0) backi += Nxh*Ny;
                }
                else if(d==2) {
                    backi -= Nxh*Ny;
                    if(z==0) backi += Nxh*Ny*Nz;
                }
                else {
                    backi -= Nxh*Ny*Nz;
                    if(t==0) backi += Nxh*Ny*Nz*Nt;
                }
#if COMPRESSED_12
            	for(int c1=0; c1<2; ++c1)
#else
            	for(int c1=0; c1<3; ++c1)
#endif
            	for(int c2=0; c2<3; ++c2) for(int ir=0; ir<2; ++ir)
            	{
                    dest[ind][2*d+1][c1][c2][ir][v] = rawsrcp[72*i+18*d+c1*6+c2*2+ir];
		    dest[ind][2*d][c1][c2][ir][v] = rawsrc[72*backi+18*d+c1*6+c2*2+ir];
		}
            }
	}
	for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
	{
	    unpack_gauge_face_dir(tid, dest, 2*dir, 1);
	}
    }
    } /* openmp */
#if TIME_GF
    timer +=  hrtimer_get();
    //if(myRank==0) printf("QPHIX_D3_create_G_from_raw: time = %.6e sec\n ", (timer)/1.e+9);
    //fflush(stdout);
#endif
    return gf;
}

/* create gauge force from raw to qphix structure */
/* convert from anti-hermition (Ic*H) in raw to hermition (H) in QPhiX */
QPHIX_D3_Force  
*QPHIX_D3_create_F_from_raw( QPHIX_D_Real *rawsrc, QPHIX_evenodd_t evenodd )
{
    MYASSERT(evenodd==QPHIX_EVENODD||evenodd==QPHIX_EVEN||evenodd==QPHIX_ODD);
#if TIME_GF
    timer = -hrtimer_get();
#endif
    QPHIX_D3_Force *gf = (QPHIX_D3_Force *)malloc(sizeof(QPHIX_D3_Force));
    MYASSERT(gf!=0x00);
    gf->heo = allocHermit();
    MYASSERT(gf->heo!=0x00);
    gf->htmp = allocHermitHelperYZT();
    MYASSERT(gf->htmp!=0x00);
    gf->parity = evenodd;
    if(evenodd == QPHIX_EVENODD) gf->parity = QPHIX_EVEN;

    Hermit *dest = (Hermit *)gf->heo;
    fptype *rawsrcp;
    if(evenodd & QPHIX_ODD) rawsrcp = rawsrc + 72*qphix_even_sites_on_node;
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    if(gf->parity == QPHIX_EVEN)
    {
	if(evenodd & QPHIX_ODD) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
        {
                pack_momentum_face_dir(tid, rawsrcp, 2*dir+1, 1);
        }

        for(int i=tid; i<qphix_even_sites_on_node; i+=nThreads)
        {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;
            int xodd = (y + z + t) & 1;

            for(int d=0; d<4; ++d)
            {
                int backi = i;
                if(d==0) {
                    backi -= 1-xodd;
                    if(xodd==0 && x==0) backi += Nxh;
                }
                else if(d==1) {
                    backi -= Nxh;
                    if(y==0) backi += Nxh*Ny;
                }
                else if(d==2) {
                    backi -= Nxh*Ny;
                    if(z==0) backi += Nxh*Ny*Nz;
                }
                else {
                    backi -= Nxh*Ny*Nz;
                    if(t==0) backi += Nxh*Ny*Nz*Nt;
                }
		for(int k=0; k<8; ++k)
		{
		    int c1 = k/4;
		    int c2 = 2;
		    int ir = k & 1;
		    if(k<2 || k==7) { 
			c2 = 1;
			if(k==7) ir = 0;
		    }
		    else if(k==6) c1 = c2 = 0;
		    ir = 1-ir;
		    dest[ind][2*d+1][k][v] = (-1.+2.*ir)*rawsrc[72*i+18*d+c1*6+c2*2+ir];
		    if(evenodd & QPHIX_ODD)
                	dest[ind][2*d][k][v] = (-1.+2.*ir)*rawsrcp[backi*72+18*d+c1*6+c2*2+ir];
		}
            }
        }
	if(evenodd & QPHIX_ODD) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
        {
                unpack_momentum_face_dir(tid, dest, 2*dir, 0);
        }
    }
    else
    {
        for(int i=tid; i<qphix_even_sites_on_node; i+=nThreads)
        {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;
            int xodd = (y + z + t + 1) & 1;

            for(int d=0; d<4; ++d)
            {
                for(int k=0; k<8; ++k)
                {
                    int c1 = k/4;
                    int c2 = 2;
		    int ir = k & 1;
                    if(k<2 || k==7) {
			c2 = 1;
			if(k==7) ir = 0;
		    }
                    else if(k==6) c1 = c2 = 0;
		    ir = 1-ir;
                    dest[ind][2*d+1][k][v] = (-1.+2.*ir)*rawsrcp[72*i+18*d+c1*6+c2*2+ir];
                }
	    }
	}
    }
    } /* openmp */
#if TIME_GF
    timer +=  hrtimer_get();
    //if(myRank==0) printf("QPHIX_D3_create_F_from_raw: time = %.6e sec\n ", (timer)/1.e+9);
    //fflush(stdout);
#endif
    return gf;
}

void 
QPHIX_D3_extract_V_to_raw( QPHIX_D_Real *dest, QPHIX_D3_ColorVector *src, QPHIX_evenodd_t evenodd )
{
#if TIME_CG
    timer = -hrtimer_get();
#endif
    int parity0 = src->parity;
    KS *kse = (KS *)src->even;
    KS *kso = (KS *)src->odd;
    if(parity0 & QPHIX_EVEN) 
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads)     
#endif
    for(int i=0; i<qphix_even_sites_on_node; ++i)
    {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

            for(int c = 0; c < 3; c++) {
		    dest[6*i+c*2] = kse[ind][c][0][v];
		    dest[6*i+c*2+1] = kse[ind][c][1][v];
	    }
    }
    if( parity0 & QPHIX_ODD )
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads)     
#endif
    for(int i=0; i<qphix_even_sites_on_node; ++i)
    {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

            for(int c = 0; c < 3; c++) {
                    dest[6*(i+qphix_even_sites_on_node)+c*2] = kso[ind][c][0][v];
                    dest[6*(i+qphix_even_sites_on_node)+c*2+1] = kso[ind][c][1][v];
            }
    }

#if TIME_CG
    timer += hrtimer_get();
    //if(myRank==0) printf("QPHIX_D3_extract_V_to_raw: time = %.6e sec\n ", (timer)/1.e+9);
    //fflush(stdout);
#endif
}

/* copy SU3 gauge fields from qphix structure to raw */
void QPHIX_D3_extract_G_to_raw( QPHIX_D_Real *dest, QPHIX_D3_GaugeField  *src, QPHIX_evenodd_t evenodd )
{
    /* check parity */
#if TIME_GF
    timer -= hrtimer_get();
#endif
    Gauge *gf = (Gauge *)src->geo;
    MYASSERT(gf!=0x00);

    QPHIX_D_Real *destp;
    if(evenodd & QPHIX_ODD) destp = dest + 72*qphix_even_sites_on_node;

    if(src->parity == QPHIX_EVEN)
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
    	int tid = omp_get_thread_num();
#else
    	int tid = 0;
#endif
	if(evenodd & QPHIX_ODD) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
	{
	    pack_gauge_face_dir(tid, gf, 2*dir, 0);
	}

	for(int i=tid; i<qphix_even_sites_on_node; i+=nThreads)
	{
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;
            int xodd = (y + z + t) & 1;

            for(int d=0; d<4; ++d)
            {
                int backi = i;
                if(d==0) {
                    backi -= 1-xodd;
                    if(xodd==0 && x==0) backi += Nxh;
                }
                else if(d==1) {
                    backi -= Nxh;
                    if(y==0) backi += Nxh*Ny;
                }
                else if(d==2) {
                    backi -= Nxh*Ny;
                    if(z==0) backi += Nxh*Ny*Nz;
                }
                else {
                    backi -= Nxh*Ny*Nz;
                    if(t==0) backi += Nxh*Ny*Nz*Nt;
                }
#if COMPRESSED_12
		for(int c1=0; c1<2; ++c1)
#else
		for(int c1=0; c1<3; ++c1)
#endif
		for(int c2=0; c2<3; ++c2) for(int ir=0; ir<2; ++ir)
		{
		    if(evenodd & QPHIX_EVEN)
		    {
			dest[72*i+18*d+6*c1+2*c2+ir] = gf[ind][2*d+1][c1][c2][ir][v];
		    }
		    if(evenodd & QPHIX_ODD)
		    {
		    	destp[72*backi+18*d+6*c1+2*c2+ir] = gf[ind][2*d][c1][c2][ir][v];
		    }
		}
#if COMPRESSED_12
		if(evenodd & QPHIX_EVEN)
		    reconstruct_gauge_third_row(&dest[72*i+18*d], 0, 0, 0);
		if(evenodd & QPHIX_ODD)
		    reconstruct_gauge_third_row(&destp[72*backi+18*d], 0, 0, 0);
#endif		
	    }
	}

	if(evenodd & QPHIX_ODD) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
		unpack_gauge_face_dir(tid, destp, 2*dir+1, 1);
    }
    else
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
    	int tid = omp_get_thread_num();
#else
    	int tid = 0;
#endif
        if(evenodd & QPHIX_EVEN) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
	{
                pack_gauge_face_dir(tid, gf, 2*dir, 1);
	}

        for(int i=tid; i<qphix_even_sites_on_node; i+=nThreads)
        {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;
            int xodd = (y + z + t + 1) & 1;

            for(int d=0; d<4; ++d)
            {
                int backi = i;
                if(d==0) {
                    backi -= 1-xodd;
                    if(xodd==0 && i==0) backi += Nxh;
                }
                else if(d==1) {
                    backi -= Nxh;
                    if(y==0) backi += Nxh*Ny;
                }
                else if(d==2) {
                    backi -= Nxh*Ny;
                    if(z==0) backi += Nxh*Ny*Nz;
                }
                else {
                    backi -= Nxh*Ny*Nz;
                    if(t==0) backi += Nxh*Ny*Nz*Nt;
                }
#if COMPRESSED_12
                for(int c1=0; c1<2; ++c1)
#else
                for(int c1=0; c1<3; ++c1)
#endif
                for(int c2=0; c2<3; ++c2) for(int ir=0; ir<2; ++ir)
                {
                    if(evenodd & QPHIX_ODD)
		    {
			destp[72*i+18*d+6*c1+2*c2+ir] = gf[ind][2*d+1][c1][c2][ir][v];
		    }
		    if(evenodd & QPHIX_EVEN)
		    {
                    	dest[72*backi+18*d+6*c1+2*c2+ir] = gf[ind][2*d][c1][c2][ir][v];
		    }
		}
#if COMPRESSED_12
		if(evenodd & QPHIX_ODD)
		    reconstruct_gauge_third_row(&destp[72*i+18*d], 0, 0, 0);
		if(evenodd & QPHIX_EVEN)
		    reconstruct_gauge_third_row(&dest[72*backi+18*d], 0, 0, 0);
#endif
	    }
	}

	if(evenodd & QPHIX_EVEN) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
		unpack_gauge_face_dir(tid, dest, 2*dir+1, 0);
    }

#if TIME_GF
    timer += hrtimer_get();
    //if(myRank==0) printf("QPHIX_D3_extract_G_to_raw: time = %.6e sec\n ", (timer)/1.e+9);
    //fflush(stdout);
#endif

}

/* copy gauge force from qphix structure to raw */
void QPHIX_D3_extract_F_to_raw( QPHIX_D_Real *dest, QPHIX_D3_Force  *src, QPHIX_evenodd_t evenodd )
{
#if TIME_GF
    timer -= hrtimer_get();
#endif
    Hermit *gf = (Hermit *)src->heo;
    MYASSERT(gf!=0x00);

    int nv[4] = {Nx, Ny, Nz, Nt};
    QPHIX_D_Real *destp = dest + 72*qphix_even_sites_on_node;

    if(src->parity == QPHIX_EVEN)
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        if(evenodd & QPHIX_ODD) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
        {
                pack_momentum_face_dir(tid, gf, 2*dir, 0);
	}

        for(int i=tid; i<qphix_even_sites_on_node; i+=nThreads)
        {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;
            int xodd = (y + z + t) & 1;

            for(int d=0; d<4; ++d)
            {
                int backi = i;
                if(d==0) {
                    backi -= 1-xodd;
                    if(xodd==0 && x==0) backi += Nxh;
                }
                else if(d==1) {
                    backi -= Nxh;
                    if(y==0) backi += Nxh*Ny;
                }
                else if(d==2) {
                    backi -= Nxh*Ny;
                    if(z==0) backi += Nxh*Ny*Nz;
                }
                else {
                    backi -= Nxh*Ny*Nz;
                    if(t==0) backi += Nxh*Ny*Nz*Nt;
                }
		for(int k=0; k<8; ++k)
		{
                    int c1 = k/4;
                    int c2 = 2;
                    int ir = k & 1;
		    if(k<2 || k==7) {
			c2 = 1;
			if(k==7) ir = 0;
		    }
                    else if(k==6) c1 = c2 = 0;
		    ir = 1-ir;
		    if(evenodd & QPHIX_EVEN) 
		    {
			dest[72*i+18*d+6*c1+2*c2+ir] = (-1.+2.*ir)*gf[ind][2*d+1][k][v];
			dest[72*i+18*d+6*c2+2*c1+ir] = gf[ind][2*d+1][k][v];
		    }
		    if(evenodd & QPHIX_ODD) 
		    {
			destp[72*backi+18*d+6*c1+2*c2+ir] = (-1.+2.*ir)*gf[ind][2*d][k][v];
			destp[72*backi+18*d+6*c2+2*c1+ir] = gf[ind][2*d][k][v];
		    }
		}
                if(evenodd & QPHIX_EVEN)
                {
                    for(int k=0; k<3; ++k)
                        dest[72*i+18*d+8*k] = 0.0;
                    dest[72*i+18*d+17] = -dest[72*i+18*d+1]-dest[72*i+18*d+9];
                }
                if(evenodd & QPHIX_ODD)
                {
                    for(int k=0; k<3; ++k)
                        destp[72*backi+18*d+8*k] = 0.0;
                    destp[72*backi+18*d+17] = -destp[72*backi+18*d+1]-destp[72*backi+18*d+9];
                }
	    }
	}

	if(evenodd & QPHIX_ODD) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
		unpack_momentum_face_dir(tid, destp, 2*dir+1, 1);
    }
    else
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        if(evenodd & QPHIX_EVEN) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
        {
                pack_momentum_face_dir(tid, gf, 2*dir, 1);
	}

        for(int i=tid; i<qphix_even_sites_on_node; i+=nThreads)
        {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;
            int xodd = (y + z + t + 1) & 1;

            for(int d=0; d<4; ++d)
            {
                int backi = i;
                if(d==0) {
                    backi -= 1-xodd;
                    if(xodd==0 && x==0) backi += Nxh;
                }
                else if(d==1) {
                    backi -= Nxh;
                    if(y==0) backi += Nxh*Ny;
                }
                else if(d==2) {
                    backi -= Nxh*Ny;
                    if(z==0) backi += Nxh*Ny*Nz;
                }
                else {
                    backi -= Nxh*Ny*Nz;
                    if(t==0) backi += Nxh*Ny*Nz*Nt;
                }
		for(int k=0; k<8; ++k)
		{
                    int c1 = k/4;
                    int c2 = 2;
		    int ir = k & 1;
                    if(k<2 || k==7) {
			c2 = 1;
			if(k==7) ir = 0;
		    }
                    else if(k==6) c1 = c2 = 0;
		    ir = 1-ir;
                    if(evenodd & QPHIX_ODD) 
		    {
			destp[72*i+18*d+6*c1+2*c2+ir] = (-1.+2.*ir)*gf[ind][2*d+1][k][v];
			destp[72*i+18*d+6*c2+2*c1+ir] = gf[ind][2*d+1][k][v];
		    }
                    if(evenodd & QPHIX_EVEN) 
		    {
			dest[72*backi+18*d+6*c1+2*c2+ir] = (-1.+2.*ir)*gf[ind][2*d][k][v];
			dest[72*backi+18*d+6*c2+2*c1+ir] = gf[ind][2*d][k][v];
		    }
		}
                if(evenodd & QPHIX_ODD)
                {
                    for(int k=0; k<3; ++k)
                        destp[72*i+18*d+8*k] = 0.0;
                    destp[72*i+18*d+17] = -destp[72*i+18*d+1]-destp[72*i+18*d+9];
                }
                if(evenodd & QPHIX_EVEN)
                {
                    for(int k=0; k<3; ++k)
                        dest[72*backi+18*d+8*k] = 0.0;
                    dest[72*backi+18*d+17] = -dest[72*backi+18*d+1]-dest[72*backi+18*d+9];
                }
	    }
	}

	if(evenodd & QPHIX_EVEN) for(int dir=0; dir<4; ++dir) if(!local_dir[dir])
		unpack_momentum_face_dir(tid, dest, 2*dir+1, 0);
    }

#if TIME_GF
    timer += hrtimer_get();
    //if(myRank==0) printf("QPHIX_D3_extract_F_to_raw: time = %.6e sec\n ", (timer)/1.e+9);
    //fflush(stdout);
#endif

}

/*
QPHIX_D3_ColorVector *
QPHIX_D3_convert_V_from_raw( QPHIX_D_Real *src, QPHIX_evenodd_t evenodd )
{
    return QPHIX_D3_create_V_from_raw( src, evenodd );
}

QPHIX_D_Real *
QPHIX_D3_convert_V_to_raw( QPHIX_D3_ColorVector *src )
{
    int nst;
    if(src->parity==QPHIX_EVENODD) nst = n_sites_on_node;
    else nst = qphix_even_sites_on_node;
    QPHIX_D_Real *R = (QPHIX_D_Real *)_mm_malloc(sizeof(QPHIX_D_Real)*nst);
    QPHIX_D3_extract_V_to_raw( R, src );
    return R;
}
*/

QPHIX_D3_FermionLinksAsqtad *
QPHIX_D3_asqtad_create_L_from_raw( QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_D_Real *fatback, QPHIX_D_Real *lngback, QPHIX_evenodd_t evenodd )
{
    string myname="QPHIX_D3_asqtad_create_L_from_raw";
#if TIME_CG
    timer = -hrtimer_get();
#endif
    QPHIX_D3_FermionLinksAsqtad *L = (QPHIX_D3_FermionLinksAsqtad *)malloc(sizeof(QPHIX_D3_FermionLinksAsqtad));
    MYASSERT(L!=0x00);
    L->lngeven = allocGauge();//initBuf(_mm_malloc((Pxyz*Vt)*sizeof(Gauge), 64), (Pxyz*Vt)*sizeof(Gauge));
    MYASSERT(L->lngeven!=0x00);
    L->lngodd = allocGauge();//initBuf(_mm_malloc((Pxyz*Vt)*sizeof(Gauge), 64), (Pxyz*Vt)*sizeof(Gauge));
    MYASSERT(L->lngodd!=0x00);
    L->fateven = allocGauge18(); //initBuf(_mm_malloc((Pxyz*Vt)*sizeof(Gauge18), 64), (Pxyz*Vt)*sizeof(Gauge18));
    MYASSERT(L->fateven!=0x00);
    L->fatodd = allocGauge18(); //initBuf(_mm_malloc((Pxyz*Vt)*sizeof(Gauge18), 64), (Pxyz*Vt)*sizeof(Gauge18));
    MYASSERT(L->fatodd!=0x00);
    Gauge *glle = (Gauge *)L->lngeven;
    Gauge *gllo = (Gauge *)L->lngodd;
    Gauge18 *gfle = (Gauge18 *)L->fateven;
    Gauge18 *gflo = (Gauge18 *)L->fatodd;

    if(evenodd & QPHIX_EVEN)
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads)     
#endif
    for(int i=0; i<qphix_even_sites_on_node; ++i)
    {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

	    for( int dir=0; dir<4; ++dir )
	    {
		int xdirf = 2*dir+1;
		int xdirb = 2*dir;
	    	for( int r=0; r<3; ++r ) for( int c=0; c<3; ++c ) 
	    	{
			gfle[ind][xdirb][r][c][0][v] = fatback[72*i+18*dir+r*6+c*2];
			gfle[ind][xdirb][r][c][1][v] = fatback[72*i+18*dir+r*6+c*2+1];
			gfle[ind][xdirf][r][c][0][v] = fat[72*i+18*dir+r*6+c*2];
			gfle[ind][xdirf][r][c][1][v] = fat[72*i+18*dir+r*6+c*2+1];
#ifdef COMPRESSED_12
		    	if(r < 2)
#endif
		    	{
			    glle[ind][xdirb][r][c][0][v] = lngback[72*i+18*dir+r*6+c*2];
			    glle[ind][xdirb][r][c][1][v] = lngback[72*i+18*dir+r*6+c*2+1];
			    glle[ind][xdirf][r][c][0][v] = lng[72*i+18*dir+r*6+c*2];
		    	    glle[ind][xdirf][r][c][1][v] = lng[72*i+18*dir+r*6+c*2+1];
		    	}
	    	}
	    }
    }
    if(evenodd & QPHIX_ODD) 
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads)     
#endif
    for(int i=0; i<qphix_even_sites_on_node; ++i)
    {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

            for( int dir=0; dir<4; ++dir )
            {
                int xdirf = 2*dir+1;
                int xdirb = 2*dir;
                for( int r=0; r<3; ++r ) for( int c=0; c<3; ++c )
                {
                        gflo[ind][xdirb][r][c][0][v] = fatback[72*(i+qphix_even_sites_on_node)+18*dir+r*6+c*2];
                        gflo[ind][xdirb][r][c][1][v] = fatback[72*(i+qphix_even_sites_on_node)+18*dir+r*6+c*2+1];
                        gflo[ind][xdirf][r][c][0][v] = fat[72*(i+qphix_even_sites_on_node)+18*dir+r*6+c*2];
                        gflo[ind][xdirf][r][c][1][v] = fat[72*(i+qphix_even_sites_on_node)+18*dir+r*6+c*2+1];
#ifdef COMPRESSED_12
                        if(r < 2)
#endif
                        {
                            gllo[ind][xdirb][r][c][0][v] = lngback[72*(i+qphix_even_sites_on_node)+18*dir+r*6+c*2];
                            gllo[ind][xdirb][r][c][1][v] = lngback[72*(i+qphix_even_sites_on_node)+18*dir+r*6+c*2+1];
                            gllo[ind][xdirf][r][c][0][v] = lng[72*(i+qphix_even_sites_on_node)+18*dir+r*6+c*2];
                            gllo[ind][xdirf][r][c][1][v] = lng[72*(i+qphix_even_sites_on_node)+18*dir+r*6+c*2+1];
                        }
                }
            }
    }

#if TIME_CG
    timer += hrtimer_get();
    //if(myRank==0) cout << myname << ": time = " << timer/1.e+9 << "sec\n" << endl;
#endif
    return L;
}

void 
reconstruct_gauge_third_row( QPHIX_D_Real *g, int t, int x, int y )
{
    for(int d=0; d<4; ++d) for(int c=0; c<3; ++c) 
    {
	g[18*d+12+c*2] = g[18*d+(c+1)%3*2] * g[18*d+6+(c+2)%3*2] - g[18*d+(c+1)%3*2+1] * g[18*d+6+(c+2)%3*2+1] - g[18*d+(c+2)%3*2] * g[18*d+6+(c+1)%3*2] + g[18*d+(c+2)%3*2+1] * g[18*d+6+(c+1)%3*2+1];
	g[18*d+12+c*2+1] = - g[18*d+(c+1)%3*2] * g[18*d+6+(c+2)%3*2+1] - g[18*d+(c+1)%3*2+1] * g[18*d+6+(c+2)%3*2] + g[18*d+(c+2)%3*2] * g[18*d+6+(c+1)%3*2+1] + g[18*d+(c+2)%3*2+1] * g[18*d+6+(c+1)%3*2];
	if( d==0 && t&1 || d==1 && (t+x)&1 || d==2 && (t+x+y)&1 ) 
	{
	    g[18*d+12+c*2] *= -1.0; g[18*d+12+c*2+1] *= -1.0; 
	}
    }
}

void 
QPHIX_D3_asqtad_extract_L_to_raw( QPHIX_D_Real *fat, QPHIX_D_Real *lng, QPHIX_D3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd )
{
    string myname = "QPHIX_D3_asqtad_extract_L_to_raw";
#if TIME_CG
    timer = -hrtimer_get();
#endif
    Gauge *glle = (Gauge *)src->lngeven;
    Gauge *gllo = (Gauge *)src->lngodd;
    Gauge18 *gfle = (Gauge18 *)src->fateven;
    Gauge18 *gflo = (Gauge18 *)src->fatodd;
    QPHIX_D_Real *fatp, *lngp;
    if(evenodd & QPHIX_ODD)
    {
	fatp = fat + qphix_even_sites_on_node*72;
	lngp = lng + qphix_even_sites_on_node*72;
    }
    if(evenodd & QPHIX_EVEN)
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads)     
#endif
    for(int i=0; i<qphix_even_sites_on_node; ++i)
    {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

            for( int dir=0; dir<4; ++dir )
	    {
		int xdir = 2*dir+1;
                for( int r=0; r<3; ++r ) for( int c=0; c<3; ++c )
                {
		    	fat[72*i+18*dir+r*6+c*2] = gfle[ind][xdir][r][c][0][v];
		    	fat[72*i+18*dir+r*6+c*2+1] = gfle[ind][xdir][r][c][1][v];
#ifdef COMPRESSED_12
		    	if( r<2 ) 
#endif
		    	{
		    	    lng[72*i+18*dir+r*6+c*2] = glle[ind][xdir][r][c][0][v];
		    	    lng[72*i+18*dir+r*6+c*2+1] = glle[ind][xdir][r][c][1][v];
		    	}
	    	}
#ifdef COMPRESSED_12
		reconstruct_gauge_third_row(lng+72*i, t, (y+z+t)&1, y);
#endif
	    }
    }
    if(evenodd & QPHIX_ODD)
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads)     
#endif
    for(int i=0; i<qphix_even_sites_on_node; ++i)
    {
            int y = i / Nxh;
            int x = i - y * Nxh;
            int z = y / Ny;
            y -= z * Ny;
            int t = z / Nz;
            z -= t * Nz;

            int x1, x2, y1, y2, z1, z2, t1, t2;
            x1 = x / Vxh; x2 = x % Vxh;
            y1 = y / Vy; y2 = y % Vy;
            z1 = z / Vz; z2 = z % Vz;
            t1 = t / Vt; t2 = t % Vt;

            int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
            int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

            for( int dir=0; dir<4; ++dir )
            {
                int xdir = 2*dir+1;
                for( int r=0; r<3; ++r ) for( int c=0; c<3; ++c )
                {
                        fatp[72*i+18*dir+r*6+c*2] = gflo[ind][xdir][r][c][0][v];
                        fatp[72*i+18*dir+r*6+c*2+1] = gflo[ind][xdir][r][c][1][v];
#ifdef COMPRESSED_12
                        if( r<2 )
#endif
                        {
                            lngp[72*i+18*dir+r*6+c*2] = gllo[ind][xdir][r][c][0][v];
                            lngp[72*i+18*dir+r*6+c*2+1] = gllo[ind][xdir][r][c][1][v];
                        }
                }
#ifdef COMPRESSED_12
		reconstruct_gauge_third_row(lngp+72*i, t, (y+z+t)^1, y);
#endif
            }
    }

#if TIME_CG
    timer += hrtimer_get();
    //if(myRank==0) cout << myname << ": time = " << timer/1.e+9 << "sec\n" << endl;
#endif
}

/*
QPHIX_D3_FermionLinksAsqtad *
QPHIX_D3_asqtad_convert_L_from_raw( QPHIX_D_Real *fat[], QPHIX_D_Real *lng[], QPHIX_evenodd_t evenodd )
{

}

void 
QPHIX_D3_asqtad_convert_L_to_raw( QPHIX_D_Real *fat[], QPHIX_D_Real *lng[], QPHIX_D3_FermionLinksAsqtad *src, QPHIX_evenodd_t evenodd )
{

}

void 
QPHIX_D3_asqtad_load_L_from_raw( QPHIX_D3_FermionLinksAsqtad *dest, QPHIX_D_Real *fat[], QPHIX_D_Real *lng[], QPHIX_evenodd_t evenodd )
{

}
*/

void 
QPHIX_D3_destroy_V( QPHIX_D3_ColorVector *V )
{
  if(V==0x00) return;
  if(V->even!=0x00) _mm_free(V->even);
  if(V->odd!=0x00) _mm_free(V->odd);
  free(V);
}

void 
QPHIX_D3_asqtad_destroy_L( QPHIX_D3_FermionLinksAsqtad *L )
{
  if(L==0x00) return;
  if(L->lngeven!=0x00) _mm_free(L->lngeven);
  if(L->lngodd!=0x00) _mm_free(L->lngodd);
  if(L->fateven!=0x00) _mm_free(L->fateven);
  if(L->fatodd!=0x00) _mm_free(L->fatodd);
  free(L);
}

void QPHIX_D3_destroy_G(QPHIX_D3_GaugeField *field)
{
    if(field==0x00) return;
    if(field->geo!=0x00) _mm_free(field->geo);
    free(field);
}

void QPHIX_D3_destroy_F(QPHIX_D3_Force *field)
{
    if(field==0x00) return;
    if(field->heo!=0x00) _mm_free(field->heo);
    if(field->htmp!=0x00) _mm_free(field->htmp);
    free(field);
    field=0x00;
}

void 
QPHIX_D3_asqtad_dslash( QPHIX_D3_FermionLinksAsqtad *asqtad,
                        QPHIX_D3_ColorVector *out,
                        QPHIX_D3_ColorVector *in,
                        QPHIX_evenodd_t parity )
{
    MYASSERT(parity != QPHIX_EVENODD);
    MYASSERT(parity & out->parity && out->parity == QPHIX_OPP_PAR(in->parity));
    if(parity==QPHIX_EVEN)
        qphix_ks_dslash( (void *)in->odd, asqtad->lngeven, asqtad->fateven, (void *)out->even, 0 );
    else
	qphix_ks_dslash( (void *)in->even, asqtad->lngodd, asqtad->fatodd, (void *)out->odd, 1 );
}

/* Convertion between SU3_Matrix & QPHIX_ColorMatrix */
void QPHIX_D3_layout_from_su3m(QPHIX_D3_ColorMatrix *dest, SU3_D_Matrix *src)
{
  int x, y, z, t;
  int x1, y1, z1, t1;
  int x2, y2, z2, t2;
  int ind, v;
  int i, c1, c2;

  QSU3M_double *dest_e = (QSU3M_double *)dest->even;
  QSU3M_double *dest_o = (QSU3M_double *)dest->odd;

#pragma omp parallel for private(x, y, z, t, x1, y1, z1, t1, x2, y2, z2, t2, ind, v)
  for(i = 0; i < qphix_even_sites_on_node; ++i)
  {
    y = i / Nxh;
    x = i - y * Nxh;
    z = y / Ny;
    y = y - z * Ny;
    t = z / Nz;
    z = z - t * Nz;

    x1 = x / Vxh;
    x2 = x % Vxh;
    y1 = y / Vy;
    y2 = y % Vy;
    z1 = z / Vz;
    z2 = z % Vz;
    t1 = t / Vt;
    t2 = t % Vt;

    ind = t2 * Pxyz + z2 * Pxy + y2 * Vxh + x2;
    v = ((t1 * nGZ + z1) * nGY + y1) * nGX + x1;

    for(c1 = 0; c1 < 3; ++c1)
    {
      for(c2 = 0; c2 < 3; ++c2)
      {
        dest_e[ind][c1][c2][0][v] = src[i].e[c1][c2].real;
        dest_e[ind][c1][c2][1][v] = src[i].e[c1][c2].imag;
        dest_o[ind][c1][c2][0][v] = src[i + qphix_even_sites_on_node].e[c1][c2].real;
        dest_o[ind][c1][c2][1][v] = src[i + qphix_even_sites_on_node].e[c1][c2].imag;
      }
    }
  }
}

void QPHIX_D3_layout_from_4su3m(QPHIX_D3_ColorMatrix *dest[], SU3_D_Matrix *src)
{
  int x, y, z, t;
  int x1, y1, z1, t1;
  int x2, y2, z2, t2;
  int ind, v;
  int i, c1, c2;

  QSU3M_double *dest_e[4];
  QSU3M_double *dest_o[4];
  
  for(int dir = 0; dir < 4; dir++){
    dest_e[dir] = (QSU3M_double *)dest[dir]->even;
    dest_o[dir] = (QSU3M_double *)dest[dir]->odd;
  }

#pragma omp parallel for private(x, y, z, t, x1, y1, z1, t1, x2, y2, z2, t2, ind, v)
  for(i = 0; i < qphix_even_sites_on_node; ++i)
  {
    y = i / Nxh;
    x = i - y * Nxh;
    z = y / Ny;
    y = y - z * Ny;
    t = z / Nz;
    z = z - t * Nz;

    x1 = x / Vxh;
    x2 = x % Vxh;
    y1 = y / Vy;
    y2 = y % Vy;
    z1 = z / Vz;
    z2 = z % Vz;
    t1 = t / Vt;
    t2 = t % Vt;

    ind = t2 * Pxyz + z2 * Pxy + y2 * Vxh + x2;
    v = ((t1 * nGZ + z1) * nGY + y1) * nGX + x1;

    for(c1 = 0; c1 < 3; ++c1)
    {
      for(c2 = 0; c2 < 3; ++c2)
      {
	for(int dir = 0; dir < 4; dir++){
	  dest_e[dir][ind][c1][c2][0][v] = src[4*i+dir].e[c1][c2].real;
	  dest_e[dir][ind][c1][c2][1][v] = src[4*i+dir].e[c1][c2].imag;
	  dest_o[dir][ind][c1][c2][0][v] =
	     src[4*(i + qphix_even_sites_on_node) + dir].e[c1][c2].real;
	  dest_o[dir][ind][c1][c2][1][v] =
	     src[4*(i + qphix_even_sites_on_node) + dir].e[c1][c2].imag;
	}
      }
    }
  }
}

void QPHIX_D3_layout_to_su3m(SU3_D_Matrix *dest, QPHIX_D3_ColorMatrix *src)
{
  int x, y, z, t;
  int x1, y1, z1, t1;
  int x2, y2, z2, t2;
  int ind, v;
  int i, c1, c2;
  
  QSU3M_double *qme = (QSU3M_double *)src->even;
  QSU3M_double *qmo = (QSU3M_double *)src->odd;
  
#pragma omp parallel for private(x, y, z, t, x1, y1, z1, t1, x2, y2, z2, t2, ind, v)
  for(i = 0; i < qphix_even_sites_on_node; i++)
  {
    y = i / Nxh;
    x = i - y * Nxh;
    z = y / Ny;
    y = y - z * Ny;
    t = z / Nz;
    z = z - t * Nz;
  
    x1 = x / Vxh;
    x2 = x % Vxh;
    y1 = y / Vy;
    y2 = y % Vy;
    z1 = z / Vz;
    z2 = z % Vz;
    t1 = t / Vt;
    t2 = t % Vt;
  
    ind = t2 * Pxyz + z2 * Pxy + y2 * Vxh + x2;
    v = ((t1 * nGZ + z1) * nGY + y1) * nGX + x1;
  
    for(c1 = 0; c1 < 3; ++c1)
    {
      for(c2 = 0; c2 < 3; ++c2)
      {
        dest[i].e[c1][c2].real = qme[ind][c1][c2][0][v];
        dest[i].e[c1][c2].imag = qme[ind][c1][c2][1][v];
        dest[i + qphix_even_sites_on_node].e[c1][c2].real = qmo[ind][c1][c2][0][v];
        dest[i + qphix_even_sites_on_node].e[c1][c2].imag = qmo[ind][c1][c2][1][v];
      }
    }
  }
}

void QPHIX_D3_layout_to_4su3m(SU3_D_Matrix *dest, QPHIX_D3_ColorMatrix *src[])
{
  int x, y, z, t;
  int x1, y1, z1, t1;
  int x2, y2, z2, t2;
  int ind, v;
  int i, c1, c2;
  
  QSU3M_double *qme[4];
  QSU3M_double *qmo[4];
  
  for(int dir = 0; dir < 4; dir++){
    qme[dir] = (QSU3M_double *)src[dir]->even;
    qmo[dir] = (QSU3M_double *)src[dir]->odd;
  }

#pragma omp parallel for private(x, y, z, t, x1, y1, z1, t1, x2, y2, z2, t2, ind, v)
  for(i = 0; i < qphix_even_sites_on_node; i++)
  {
    y = i / Nxh;
    x = i - y * Nxh;
    z = y / Ny;
    y = y - z * Ny;
    t = z / Nz;
    z = z - t * Nz;
  
    x1 = x / Vxh;
    x2 = x % Vxh;
    y1 = y / Vy;
    y2 = y % Vy;
    z1 = z / Vz;
    z2 = z % Vz;
    t1 = t / Vt;
    t2 = t % Vt;
  
    ind = t2 * Pxyz + z2 * Pxy + y2 * Vxh + x2;
    v = ((t1 * nGZ + z1) * nGY + y1) * nGX + x1;
  
    for(c1 = 0; c1 < 3; ++c1)
    {
      for(c2 = 0; c2 < 3; ++c2)
      {
	for(int dir = 0; dir < 4; dir++){
	  dest[4*i + dir].e[c1][c2].real = qme[dir][ind][c1][c2][0][v];
	  dest[4*i + dir].e[c1][c2].imag = qme[dir][ind][c1][c2][1][v];
	  dest[4*(i + qphix_even_sites_on_node) + dir].e[c1][c2].real =
	    qmo[dir][ind][c1][c2][0][v];
	  dest[4*(i + qphix_even_sites_on_node) + dir].e[c1][c2].imag =
	    qmo[dir][ind][c1][c2][1][v];
	}
      }
    }
  }
}

/* Creation of HISQ fermion links */
QPHIX_D3_FermionLinksHisq *
QPHIX_D3_hisq_create_L_from_4su3m(void *U_links, void *V_links, void *W_unitlinks, void *Y_unitlinks,
				  QPHIX_evenodd_t parity)
{
  QPHIX_D3_FermionLinksHisq *qphix_links = 
    (QPHIX_D3_FermionLinksHisq *)malloc(sizeof(QPHIX_D3_FermionLinksHisq));
  MYASSERT( qphix_links != NULL );

  for(int i=0; i<4; i++) {
    qphix_links->U_links[i] = QPHIX_D3_create_M(parity);
    qphix_links->V_links[i] = QPHIX_D3_create_M(parity);
    qphix_links->W_unitlinks[i] = QPHIX_D3_create_M(parity);
    qphix_links->Y_unitlinks[i] = QPHIX_D3_create_M(parity);
  }

  QPHIX_D3_layout_from_4su3m(qphix_links->U_links, (SU3_D_Matrix*)U_links);
  QPHIX_D3_layout_from_4su3m(qphix_links->V_links, (SU3_D_Matrix*)V_links);
  QPHIX_D3_layout_from_4su3m(qphix_links->W_unitlinks, (SU3_D_Matrix*)W_unitlinks);
  QPHIX_D3_layout_from_4su3m(qphix_links->Y_unitlinks, (SU3_D_Matrix*)Y_unitlinks);

  return qphix_links;
}


static void
destroy_4M(QPHIX_D3_ColorMatrix **M){
  if(M == NULL)return;
  for(int dir = 0; dir < 4; dir++)
    if(M[dir] != NULL)QPHIX_D3_destroy_M( M[dir] );
}

void 
QPHIX_D3_hisq_destroy_L( QPHIX_D3_FermionLinksHisq *L )
{
    if(L==0x00) return;
    if(L->U_links!=0x00)destroy_4M(L->U_links);
    if(L->V_links!=0x00)destroy_4M(L->V_links);
    if(L->W_unitlinks!=0x00)destroy_4M(L->W_unitlinks);
    if(L->Y_unitlinks!=0x00)destroy_4M(L->Y_unitlinks);
    free(L);
}

/* Convert SU3_AntiHermitMat to QPHIX_ColorMatrix */
//void QPHIX_D3_layout_from_anti_hermitmat(QPHIX_D3_ColorMatrix *dest[], void *ptr, int SZ)
QPHIX_D3_ColorMatrix **
QPHIX_D3_create_F_from_anti_hermitmat( void *ptr, int SZ )
{
  int x, y, z, t;
  int x1, y1, z1, t1;
  int x2, y2, z2, t2;
  int ind, v;
  int i, c1, c2;

  QPHIX_D3_ColorMatrix **dest = 
    (QPHIX_D3_ColorMatrix **)malloc(sizeof(QPHIX_D3_ColorMatrix *)*4);
  MYASSERT( dest != NULL );

  for(int dir = 0; dir < 4; dir++){
    dest[dir] = QPHIX_D3_create_M(QPHIX_EVENODD);
    MYASSERT( dest[dir] != NULL );
  }
  
  QSU3M_double *dest_e[4];
  QSU3M_double *dest_o[4];
  
  for(int dir = 0; dir < 4; dir++){
    dest_e[dir] = (QSU3M_double *)dest[dir]->even;
    dest_o[dir] = (QSU3M_double *)dest[dir]->odd;
  }
  
#pragma omp parallel for private(x, y, z, t, x1, y1, z1, t1, x2, y2, z2, t2, ind, v)
  for(i = 0; i < qphix_even_sites_on_node; i++)
  {
    y = i / Nxh;
    x = i - y * Nxh;
    z = y / Ny;
    y = y - z * Ny;
    t = z / Nz;
    z = z - t * Nz;
    
    x1 = x / Vxh;
    x2 = x % Vxh;
    y1 = y / Vy;
    y2 = y % Vy;
    z1 = z / Vz;
    z2 = z % Vz;
    t1 = t / Vt;
    t2 = t % Vt;
    
    ind = t2 * Pxyz + z2 * Pxy + y2 * Vxh + x2;
    v = ((t1 * nGZ + z1) * nGY + y1) * nGX + x1;
    SU3_D_AntiHermitMat *ahe = (SU3_D_AntiHermitMat*)(ptr +  i*SZ);
    SU3_D_AntiHermitMat *aho = (SU3_D_AntiHermitMat*)(ptr + (i + qphix_even_sites_on_node)*SZ);
    
   for(int dir = 0; dir < 4; dir++){
    dest_e[dir][ind][0][0][0][v] =  0.0;
    dest_e[dir][ind][0][0][1][v] =  ahe[dir].m00im;
    dest_e[dir][ind][1][1][0][v] =  0.0;
    dest_e[dir][ind][1][1][1][v] =  ahe[dir].m11im;
    dest_e[dir][ind][2][2][0][v] =  0.0;
    dest_e[dir][ind][2][2][1][v] =  ahe[dir].m22im;
    
    dest_e[dir][ind][0][1][0][v] =  ahe[dir].m01.real;
    dest_e[dir][ind][0][1][1][v] =  ahe[dir].m01.imag;
    dest_e[dir][ind][1][0][0][v] = -ahe[dir].m01.real;
    dest_e[dir][ind][1][0][1][v] =  ahe[dir].m01.imag;
    dest_e[dir][ind][0][2][0][v] =  ahe[dir].m02.real;
    dest_e[dir][ind][0][2][1][v] =  ahe[dir].m02.imag;
    dest_e[dir][ind][2][0][0][v] = -ahe[dir].m02.real;
    dest_e[dir][ind][2][0][1][v] =  ahe[dir].m02.imag;
    dest_e[dir][ind][1][2][0][v] =  ahe[dir].m12.real;
    dest_e[dir][ind][1][2][1][v] =  ahe[dir].m12.imag;
    dest_e[dir][ind][2][1][0][v] = -ahe[dir].m12.real;
    dest_e[dir][ind][2][1][1][v] =  ahe[dir].m12.imag;
    
    dest_o[dir][ind][0][0][0][v] =  0.0;
    dest_o[dir][ind][0][0][1][v] =  aho[dir].m00im;
    dest_o[dir][ind][1][1][0][v] =  0.0;
    dest_o[dir][ind][1][1][1][v] =  aho[dir].m11im;
    dest_o[dir][ind][2][2][0][v] =  0.0;
    dest_o[dir][ind][2][2][1][v] =  aho[dir].m22im;
    
    dest_o[dir][ind][0][1][0][v] =  aho[dir].m01.real;
    dest_o[dir][ind][0][1][1][v] =  aho[dir].m01.imag;
    dest_o[dir][ind][1][0][0][v] = -aho[dir].m01.real;
    dest_o[dir][ind][1][0][1][v] =  aho[dir].m01.imag;
    dest_o[dir][ind][0][2][0][v] =  aho[dir].m02.real;
    dest_o[dir][ind][0][2][1][v] =  aho[dir].m02.imag;
    dest_o[dir][ind][2][0][0][v] = -aho[dir].m02.real;
    dest_o[dir][ind][2][0][1][v] =  aho[dir].m02.imag;
    dest_o[dir][ind][1][2][0][v] =  aho[dir].m12.real;
    dest_o[dir][ind][1][2][1][v] =  aho[dir].m12.imag;
    dest_o[dir][ind][2][1][0][v] = -aho[dir].m12.real;
    dest_o[dir][ind][2][1][1][v] =  aho[dir].m12.imag;
   }
  }
  return dest;
}

/* Convert QPHIX_ColorMatrix to SU3_AntiHermitMat */
void QPHIX_D3_layout_to_anti_hermitmat(void *ptr, QPHIX_D3_ColorMatrix *src[], int SZ)
{
  int x, y, z, t;
  int x1, y1, z1, t1;
  int x2, y2, z2, t2;
  int ind, v;
  int i, c1, c2;
  
  QSU3M_double *src_e[4];
  QSU3M_double *src_o[4];
  
  for(int dir = 0; dir < 4; dir++){
    src_e[dir] = (QSU3M_double *)src[dir]->even;
    src_o[dir] = (QSU3M_double *)src[dir]->odd;
  }
  
#pragma omp parallel for private(x, y, z, t, x1, y1, z1, t1, x2, y2, z2, t2, ind, v)
  for(i = 0; i < qphix_even_sites_on_node; ++i)
  {
    y = i / Nxh;
    x = i - y * Nxh;
    z = y / Ny;
    y = y - z * Ny;
    t = z / Nz;
    z = z - t * Nz;
    
    x1 = x / Vxh;
    x2 = x % Vxh;
    y1 = y / Vy;
    y2 = y % Vy;
    z1 = z / Vz;
    z2 = z % Vz;
    t1 = t / Vt;
    t2 = t % Vt;
    
    ind = t2 * Pxyz + z2 * Pxy + y2 * Vxh + x2;
    v = ((t1 * nGZ + z1) * nGY + y1) * nGX + x1;
    
    SU3_D_AntiHermitMat *ahe = (SU3_D_AntiHermitMat*)(ptr +  i*SZ);
    SU3_D_AntiHermitMat *aho = (SU3_D_AntiHermitMat*)(ptr + (i + qphix_even_sites_on_node)*SZ);
    
   for(int dir = 0; dir < 4; dir++){
    ahe[dir].m00im    = src_e[dir][ind][0][0][1][v];
    ahe[dir].m11im    = src_e[dir][ind][1][1][1][v];
    ahe[dir].m22im    = src_e[dir][ind][2][2][1][v];
    ahe[dir].m01.imag = src_e[dir][ind][0][1][1][v];
    ahe[dir].m01.real = src_e[dir][ind][0][1][0][v];
    ahe[dir].m02.imag = src_e[dir][ind][0][2][1][v];
    ahe[dir].m02.real = src_e[dir][ind][0][2][0][v];
    ahe[dir].m12.imag = src_e[dir][ind][1][2][1][v];
    ahe[dir].m12.real = src_e[dir][ind][1][2][0][v];
    
    aho[dir].m00im    = src_o[dir][ind][0][0][1][v];
    aho[dir].m11im    = src_o[dir][ind][1][1][1][v];
    aho[dir].m22im    = src_o[dir][ind][2][2][1][v];
    aho[dir].m01.imag = src_o[dir][ind][0][1][1][v];
    aho[dir].m01.real = src_o[dir][ind][0][1][0][v];
    aho[dir].m02.imag = src_o[dir][ind][0][2][1][v];
    aho[dir].m02.real = src_o[dir][ind][0][2][0][v];
    aho[dir].m12.imag = src_o[dir][ind][1][2][1][v];
    aho[dir].m12.real = src_o[dir][ind][1][2][0][v];
   }
  }
}

void QPHIX_D3_reset_M(QPHIX_D3_ColorMatrix *M2, qphix_su3_matrix* M1)
{
  QSU3M_double *m2e, *m2o;
  int i, v, ind;
  m2e = (QSU3M_double *)M2->even;
  m2o = (QSU3M_double *)M2->odd;

#pragma omp parallel for private(v, ind)
  for(i=0; i<qphix_fused_sites_on_node; i++)
  {
    for(v=0; v<VECLEN; v++)
    {
      ind = (i * VECLEN * 2) + (v * 2);

      m2e[i][0][0][0][v] = M1[ind][0][0].real;
      m2e[i][0][0][1][v] = M1[ind][0][0].imag;
      m2e[i][0][1][0][v] = M1[ind][0][1].real;
      m2e[i][0][1][1][v] = M1[ind][0][1].imag;
      m2e[i][0][2][0][v] = M1[ind][0][2].real;
      m2e[i][0][2][1][v] = M1[ind][0][2].imag;
      m2e[i][1][0][0][v] = M1[ind][1][0].real;
      m2e[i][1][0][1][v] = M1[ind][1][0].imag;
      m2e[i][1][1][0][v] = M1[ind][1][1].real;
      m2e[i][1][1][1][v] = M1[ind][1][1].imag;
      m2e[i][1][2][0][v] = M1[ind][1][2].real;
      m2e[i][1][2][1][v] = M1[ind][1][2].imag;
      m2e[i][2][0][0][v] = M1[ind][2][0].real;
      m2e[i][2][0][1][v] = M1[ind][2][0].imag;
      m2e[i][2][1][0][v] = M1[ind][2][1].real;
      m2e[i][2][1][1][v] = M1[ind][2][1].imag;
      m2e[i][2][2][0][v] = M1[ind][2][2].real;
      m2e[i][2][2][1][v] = M1[ind][2][2].imag;

      ind++;

      m2o[i][0][0][0][v] = M1[ind][0][0].real;
      m2o[i][0][0][1][v] = M1[ind][0][0].imag;
      m2o[i][0][1][0][v] = M1[ind][0][1].real;
      m2o[i][0][1][1][v] = M1[ind][0][1].imag;
      m2o[i][0][2][0][v] = M1[ind][0][2].real;
      m2o[i][0][2][1][v] = M1[ind][0][2].imag;
      m2o[i][1][0][0][v] = M1[ind][1][0].real;
      m2o[i][1][0][1][v] = M1[ind][1][0].imag;
      m2o[i][1][1][0][v] = M1[ind][1][1].real;
      m2o[i][1][1][1][v] = M1[ind][1][1].imag;
      m2o[i][1][2][0][v] = M1[ind][1][2].real;
      m2o[i][1][2][1][v] = M1[ind][1][2].imag;
      m2o[i][2][0][0][v] = M1[ind][2][0].real;
      m2o[i][2][0][1][v] = M1[ind][2][0].imag;
      m2o[i][2][1][0][v] = M1[ind][2][1].real;
      m2o[i][2][1][1][v] = M1[ind][2][1].imag;
      m2o[i][2][2][0][v] = M1[ind][2][2].real;
      m2o[i][2][2][1][v] = M1[ind][2][2].imag;
    }
  }
}

//Expose an array of su3 matrix from QPHIX_ColorMatrix
qphix_su3_matrix* QPHIX_D3_expose_M(QPHIX_D3_ColorMatrix* M1)
{
  qphix_su3_matrix *M2;
  QSU3M_double *m1e, *m1o;
  int i, v, ind;
  m1e = (QSU3M_double *)M1->even;
  m1o = (QSU3M_double *)M1->odd;

  //malloc for M2
  M2 = (qphix_su3_matrix *) _mm_malloc(sizeof(qphix_su3_matrix) * qphix_even_sites_on_node * 2, VECLENBYTES);

#pragma omp parallel for private(v, ind)
  for(i=0; i<qphix_fused_sites_on_node; i++)
  {
    for(v=0; v<VECLEN; v++)
    {
      ind = (i * VECLEN * 2) + (v * 2);

      M2[ind][0][0].real = m1e[i][0][0][0][v];
      M2[ind][0][0].imag = m1e[i][0][0][1][v];
      M2[ind][0][1].real = m1e[i][0][1][0][v];
      M2[ind][0][1].imag = m1e[i][0][1][1][v];
      M2[ind][0][2].real = m1e[i][0][2][0][v];
      M2[ind][0][2].imag = m1e[i][0][2][1][v];
      M2[ind][1][0].real = m1e[i][1][0][0][v];
      M2[ind][1][0].imag = m1e[i][1][0][1][v];
      M2[ind][1][1].real = m1e[i][1][1][0][v];
      M2[ind][1][1].imag = m1e[i][1][1][1][v];
      M2[ind][1][2].real = m1e[i][1][2][0][v];
      M2[ind][1][2].imag = m1e[i][1][2][1][v];
      M2[ind][2][0].real = m1e[i][2][0][0][v];
      M2[ind][2][0].imag = m1e[i][2][0][1][v];
      M2[ind][2][1].real = m1e[i][2][1][0][v];
      M2[ind][2][1].imag = m1e[i][2][1][1][v];
      M2[ind][2][2].real = m1e[i][2][2][0][v];
      M2[ind][2][2].imag = m1e[i][2][2][1][v];

      ind++;

      M2[ind][0][0].real = m1o[i][0][0][0][v];
      M2[ind][0][0].imag = m1o[i][0][0][1][v];
      M2[ind][0][1].real = m1o[i][0][1][0][v];
      M2[ind][0][1].imag = m1o[i][0][1][1][v];
      M2[ind][0][2].real = m1o[i][0][2][0][v];
      M2[ind][0][2].imag = m1o[i][0][2][1][v];
      M2[ind][1][0].real = m1o[i][1][0][0][v];
      M2[ind][1][0].imag = m1o[i][1][0][1][v];
      M2[ind][1][1].real = m1o[i][1][1][0][v];
      M2[ind][1][1].imag = m1o[i][1][1][1][v];
      M2[ind][1][2].real = m1o[i][1][2][0][v];
      M2[ind][1][2].imag = m1o[i][1][2][1][v];
      M2[ind][2][0].real = m1o[i][2][0][0][v];
      M2[ind][2][0].imag = m1o[i][2][0][1][v];
      M2[ind][2][1].real = m1o[i][2][1][0][v];
      M2[ind][2][1].imag = m1o[i][2][1][1][v];
      M2[ind][2][2].real = m1o[i][2][2][0][v];
      M2[ind][2][2].imag = m1o[i][2][2][1][v];
    }
  }

  return M2;
}

