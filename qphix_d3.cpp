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
#if TIME_CG
    timer +=  hrtimer_get();
    fflush(stdout);
#endif
    return V;
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
#endif
}

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
#endif
}

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

