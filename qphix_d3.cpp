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
//extern int sites_on_node;

static hrtimer_t timer = 0;

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
    V=0x00;
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
    L=0x00;
}

void QPHIX_D3_destroy_G(QPHIX_D3_GaugeField *field)
{
    if(field==0x00) return;
    if(field->geo!=0x00) _mm_free(field->geo);
    free(field);
    field=0x00;
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

