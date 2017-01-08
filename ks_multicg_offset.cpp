/******* ks_multicg_offset.cpp - multi-mass CG for SU3/fermions ****/

/* Multi-mass CG inverter for staggered fermions */

/* Based on B. Jegerlehner, hep-lat/9612014.
   See also A. Frommer, S. G\"usken, T. Lippert, B. N\"ockel,"
   K. Schilling, Int. J. Mod. Phys. C6 (1995) 627. 
*/
#include <string.h>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <cassert>

#include "ks_config.h"
#include "ks_globals.h"
#include "ks_d_congrad_fn.h"
#include "ks_long_dslash.h"
#include "ks_blas_utils_c.h"
#include "ks_boundary.h"

#define DEBUG 0
#define TIME_CG 1

#include <sys/time.h>
#if TIME_CG
#include "hrtimer.h"   /* High resolution timer */
#include "ks_cg_perf_profiler.h"
#endif

#include "qphix_internal.h"

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
#elif QPHIX_PrecisionInt == 2
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

extern int qphix_even_sites_on_node;

/******************************** BLAS Kernels ********************************/
static double
calc_source_norm (fptype *restrict spinor)
{
    fptype source_norm = 0.0;
#ifdef _OPENMP
#pragma omp parallel for num_threads (nThreads) reduction(+:source_norm)    
#endif
    for(int i = 0; i < num_floats_in_ks_array; ++i) {
        source_norm += spinor[i]*spinor[i];
    }
    return double(source_norm);
}

static void
copyKS( fptype *restrict dst
      , fptype*restrict src
      )
{
#pragma omp parallel for num_threads (nThreads)
    for(int i=0; i<num_floats_in_ks_array; ++i)
    {
	dst[i] = src[i];
    }
}

static void
copyKSn( fptype **restrict pm
        , fptype*restrict src
        , int num_offsets
        )
{
#pragma omp parallel for num_threads (nThreads) collapse(2)
    for(int j=0; j<num_offsets; ++j) for(int i=0; i<num_floats_in_ks_array; ++i)
    {
        pm[j][i] = src[i];
    }
}

static double
calc_pkp ( fptype *restrict ttt
         , fptype *restrict cg_p
         , fptype msq_x4)
{
    fptype pkp = 0.0;
    //int block=(num_floats_in_ks_array/(nThreads*VECLEN))*VECLEN;
#ifdef _OPENMP
#pragma omp parallel for num_threads (nThreads) reduction(+:pkp)    
#endif
    //for(int i0 = 0; i0 < num_floats_in_ks_array; i0+=block)
    //for(int i1=i0; i1<i0+block && i1<num_floats_in_ks_array; i1+=VECLEN ) 
    //for(int i=i1; i<i1+VECLEN; ++i) 
    for(int i=0; i<num_floats_in_ks_array; ++i)
    {
        scalar_mult_add_su3_vector2(ttt[i], cg_p[i], msq_x4, ttt[i]);
        pkp += cg_p[i] * ttt[i];
    }
    return double(pkp);
}

/*
static double
calc_rsq ( fptype *restrict ttt
         , const fptype *restrict t_dest
         , const fptype *restrict t_src
         , fptype *restrict resid
         , fptype *restrict cg_p
         , fptype msq_x4)
{
    fptype rsq = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:rsq)
    for(int i = 0; i < num_floats_in_ks_array; ++i) {
        scalar_mult_add_su3_vector2(ttt[i], t_dest[i], -msq_x4, ttt[i]);
        add_su3_vector2(t_src[i], ttt[i], resid[i]);
        cg_p[i] = resid[i];
        rsq += resid[i] * resid[i];
    }
    return double(rsq);
}
*/

static inline double
calc_rsq2 ( fptype *restrict resid
          , fptype a
          , const fptype *restrict ttt
          )
{
    fptype rsq = 0.0;
    //int block=(num_floats_in_ks_array/(nThreads*VECLEN))*VECLEN;
#ifdef _OPENMP
#pragma omp parallel for num_threads (nThreads) reduction(+:rsq) 
#endif
    /* for(int j=0; j<num_offsets; ++j)*/ 
    //for(int i0 = 0; i0 < num_floats_in_ks_array; i0+=block) 
    //for( int i1=i0; i1<i0+block && i1<num_floats_in_ks_array; i1+=VECLEN )
    //for(int i=i1; i<i1+VECLEN; ++i)
//    for(int i=i0; i<i0+VECLEN; ++i)
    for(int i=0; i<num_floats_in_ks_array; ++i)
    {
        //scalar_mult_add_su3_vector2(((fptype*)t_dest)[i], ((fptype*)cg_p)[i], a[j], ((fptype*)t_dest)[i]);
        //if(j==j_low) {
          scalar_mult_add_su3_vector2(resid[i], ttt[i], a, resid[i]);
          rsq += resid[i] * resid[i];
        //}
    }
    return (double) rsq;
}

static void
update_dest( KS **restrict t_dest
	    , int num_offsets
	    , const KS **restrict cg_p
	    , fptype *a
	   ) 
{
    int block=(num_floats_in_ks_array/(nThreads*VECLEN))*VECLEN;
#ifdef _OPENMP
#pragma omp parallel for num_threads (nThreads) collapse(2)
#endif
    for(int j=0; j<num_offsets; ++j) for(int i0=0; i0<num_floats_in_ks_array; i0+=block)
    for( int i1=i0; i1<i0+block && i1<num_floats_in_ks_array; i1+=VECLEN )
    for(int i=i1; i<i1+VECLEN; ++i)
    {
	scalar_mult_add_su3_vector2(((fptype*)t_dest[j])[i], ((fptype*)cg_p[j])[i], a[j], ((fptype*)t_dest[j])[i]);
    }
}

static void
update_pm ( const fptype *restrict resid
          , KS **restrict t_dest /*, fptype *restrict ttt*/
          , KS **restrict pm
	  , fptype *restrict beta_i
          , fptype *restrict zeta_ip1
          , fptype *restrict alpha
          , int num_offsets
          )
{
    int block=(num_floats_in_ks_array/(nThreads*VECLEN))*VECLEN;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads (nThreads) collapse(2)
#endif
    for(int j=0; j<num_offsets; ++j) for(int i0=0; i0<num_floats_in_ks_array; i0+=block) //for(int j=0; j<num_offsets; ++j)
    for( int i1=i0; i1<i0+block && i1<num_floats_in_ks_array; i1+=VECLEN )
    for(int i=i1; i<i1+VECLEN; ++i) 
    {
	scalar_mult_add_su3_vector2(((fptype*)t_dest[j])[i], ((fptype*)pm[j])[i], beta_i[j], ((fptype*)t_dest[j])[i]);
	fptype ttt;
        scalar_mult_su3_vector2(resid[i], zeta_ip1[j], ttt);
        scalar_mult_add_su3_vector2(ttt, ((fptype*)pm[j])[i], alpha[j], ((fptype*)pm[j])[i]);
    }
}

#undef QPHIX_fptype_asqtad_invert_multi
#undef QPHIX_fptype_FermionLinksAsqtad
#undef QPHIX_fptype_Real
#undef QPHIX_fptype_ColorVector

#if QPHIX_PrecisionInt==1
#define QPHIX_fptype_asqtad_invert_multi QPHIX_F3_asqtad_invert_multi
#define QPHIX_fptype_FermionLinksAsqtad QPHIX_F3_FermionLinksAsqtad
#define QPHIX_fptype_Real QPHIX_F_Real
#define QPHIX_fptype_ColorVector QPHIX_F3_ColorVector
#elif QPHIX_PrecisionInt==2
#define QPHIX_fptype_asqtad_invert_multi QPHIX_D3_asqtad_invert_multi
#define QPHIX_fptype_FermionLinksAsqtad QPHIX_D3_FermionLinksAsqtad
#define QPHIX_fptype_Real QPHIX_D_Real
#define QPHIX_fptype_ColorVector QPHIX_D3_ColorVector
#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif

int 
/*
qphix_ks_multicg_offset ( void *t_src_arg
                        , void **t_dest_arg
                        , struct QuarkInvertControl *qic
                        , fptype *mass
                        , int num_offsets
                        , void *gll_arg[2]
                        , void *gfl_arg[2]
                        )
*/
QPHIX_fptype_asqtad_invert_multi( QPHIX_info_t *info,
                              QPHIX_fptype_FermionLinksAsqtad *asqtad,
			      QPHIX_invert_arg_t *inv_arg,
			      QPHIX_resid_arg_t *res_arg[],
                              QPHIX_fptype_Real *mass,
                              int num_offsets,
			      QPHIX_fptype_ColorVector *out[],
			      QPHIX_fptype_ColorVector *in )
{
    if(myRank==0) printf("QPHIX_fptype_asqtad_invert_multi\n");
    fflush(stdout);
/*
#if QPHIX_PrecisionInt==2
#warning "QPHIX asqtad inverter: Right now support single precision only" 
    return 0;
#endif
*/

  if( num_offsets==0 ) {
    for(int i=0; i<num_offsets; ++i)
    {
	res_arg[i]->final_restart = 0;
	res_arg[i]->final_iter = 0;
    }
    return(0);
  }

  char myname[25] = "qphix_ks_multicg_offset";
#if TIME_CG
    hrtimer_t t = hrtimer_get();
    CGPerfData cgPerfData;
#endif
    MYASSERT(in->parity==out[0]->parity&&in->parity==inv_arg->parity);

  /* Initialize gBar and phae */
/*
    int sz = 1, sy = 1;
    if(n_threads_per_core == 1) sz = 1;
    else if(n_threads_per_core == 2) sz = 2;
    else if(n_threads_per_core == 4) {
        sy = 2;
        sz = 2;
    }
    else {
        printf("Threads per core must be set to 1, 2 or 4. Value is %d\n", n_threads_per_core);
        return QPHIX_FAIL;
    }

	printf("Initialize gBar and phase\n");
	printf("Vxh = %d NCores = %d minCt = %d n_threads_per_core = %d\n", Vxh, NCores, minCt, n_threads_per_core);
	fflush(stdout);
                gBar = new Barrier(NCores,n_threads_per_core);
                if(gBar==0x00) return QPHIX_MEM_ERROR;
                phaser = new Phaser(Vxh, Vy, Vz, Vt, Vxh
                      , BY, BZ, NCores, sy, sz, minCt);
                if(phaser==0x00) return QPHIX_MEM_ERROR;

    nThreads = NCores * n_threads_per_core;

#ifdef _OPENMP
#pragma omp parallel default (shared) num_threads(nThreads)
    {
        int tid = omp_get_thread_num();
        gBar->init(tid);
        phaser->init(tid);
    }
#else
    if(nThreads != 1) {printf("No. of threads set to %d while not using OpenMP\n", nThreads); }
    phaser->init(0);
#endif
*/

    int parity0 = inv_arg->parity;
    Gauge *gll[2];
    Gauge18 *gfl[2];
    KS *t_src;
    if(parity0!=QPHIX_ODD)
    {
   	t_src = (KS*)in->even;
    }
    else
    {
	t_src = (KS*)in->odd;
    }

/*
    for(int i=0; i<num_floats_in_ks_array; ++i)
	printf("t_src[%d] = %f\n", i, ((fptype*)t_src)[i]);
*/
    KS *t_dest[MAXNDMASS];
    MYASSERT(MAXNDMASS>=num_offsets);
    if(parity0!=QPHIX_ODD) for(int i=0; i<num_offsets; ++i)
    	t_dest[i] = (KS*)out[i]->even;
    else for(int i=0; i<num_offsets; ++i)
	t_dest[i] = (KS*)out[i]->odd;
    gll[0] = (Gauge*)asqtad->lngeven;
    gll[1] = (Gauge*)asqtad->lngodd;

    gfl[0] = (Gauge18*)asqtad->fateven;
    gfl[1] = (Gauge18*)asqtad->fatodd;

    int iteration, total_iters;                      /* counter for iterations */
#ifdef FEWSUMS
    double actual_rsq;                  /* rsq from actual summation of resid */
    double c_tr,c_tt,tempsum[4];        /* Re<resid|ttt>, <ttt|ttt> */
#endif
    double rsq = 0., relrsq = 0.;         /* resid**2, rel resid*2 */
    double oldrsq, pkp;                 /* last resid*2, pkp = cg_p.K.cg_p */
    double source_norm;                 /* squared magnitude of source vector */
    int nrestart = 0;                       /* Restart counter */

    /* temporary storages */
    KS *ttt, *ttt_temp, *resid;
    int j, j_low;
    fptype * restrict shifts;
    fptype * restrict zeta_i, * restrict zeta_im1, * restrict zeta_ip1;
    fptype * restrict beta_i, * restrict beta_im1, * restrict alpha;
    int * restrict finished;
    KS *pm[MAXNDMASS];

    /* Unpack structure */
    int niter          = inv_arg->max;
    int max_restarts   = inv_arg->nrestart;
    fptype rsqmin      = res_arg[0]->resid;
    fptype relrsqmin   = res_arg[0]->relresid;
    int max_cg = max_restarts * niter;  /* Maximum number of iterations */
    if(myRank==0)
    {
	printf("niter = %d max_restarts = %d max_cg = %d\n", niter, max_restarts, max_cg);
	fflush(stdout);
    }

    /* Convert MILC parity to QPhiX's parity. Unlike MILC, QPhiX only has 0 or 1 
 *      * as parity. */
    int parity, otherparity;
    /* if we want both parities, we will do even first. */
    /*
    switch(parity0){
    case(QPHIX_EVEN): parity=QPHIX_EVEN; otherparity=QPHIX_ODD; break;
    case(QPHIX_ODD):  parity=QPHIX_ODD; otherparity=QPHIX_EVEN; break;
    case(QPHIX_EVENODD):  parity=QPHIX_EVEN; otherparity=QPHIX_ODD; break;
    }
    */
    parity =  parity0 & 1;
    otherparity = 1 - parity;
    /* Allocate the temps ttt, cg_p, resid */
#if TIME_CG
    hrtimer_t t_misc = hrtimer_get();
#endif
    ttt      = (KS*)allocKS();
    MYASSERT(ttt!=0x00);
    ttt_temp = (KS*)allocKS();
    MYASSERT(ttt_temp!=0x00);
    resid    = (KS*)allocKS();
    MYASSERT(resid!=0x00);
    for(j=0; j<num_offsets; ++j)
    {
	pm[j] = (KS *)allocKS();
	MYASSERT(pm[j]!=0x00);
    }

    shifts = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_i = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_im1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_ip1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    beta_i = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    beta_im1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    alpha = (fptype * restrict )malloc(num_offsets*sizeof(fptype));

    finished = (int *)malloc(sizeof(int)*num_offsets);
#if TIME_CG
    //hrtimer_t t = hrtimer_get();
    //std::cout << "allocating time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif
    double offset_low = 1.0e+20;
    j_low = -1;
    for(j=0;j<num_offsets;j++){
      shifts[j] = mass[j]*mass[j]*4.0;
      if (shifts[j] < offset_low){
        offset_low = shifts[j];
        j_low = j;
      }
    }

    for(j=0;j<num_offsets;j++){
      if( j!=j_low )shifts[j] -= shifts[j_low];
    }

  iteration = 0;
  total_iters = 0;

    struct timeval tv;
    gettimeofday( &tv, NULL );
    /* Initialize info */
    info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec;
    info->final_flop = 0.;

/*
  double g_norm = 0.0;
  for(int j=0; j<2; ++j)
  { g_norm = 0.0;
    for(int i=0; i<Pxyz*Vt*sizeof(Gauge)/sizeof(fptype); ++i)
      if((i/(sizeof(Gauge)/sizeof(fptype)/8)) & 1) g_norm += ((fptype*)gll[j])[i] * ((fptype*)gll[j])[i];
    printf("gll_norm[%d] = %f\n", j, g_norm);
  }
 
  for(int j=0; j<2; ++j) {
  g_norm = 0.0;
  for(int i=0; i<Pxyz*Vt*sizeof(Gauge)/sizeof(fptype); ++i)
      if((i/(sizeof(Gauge)/sizeof(fptype)/8)+1) & 1) g_norm += ((fptype*)gll[j])[i] * ((fptype*)gll[j])[i];
  printf("gllback_norm[%d] = %f\n", j, g_norm);
  }

  for(int j=0; j<2; ++j)
  { g_norm = 0.0;
    for(int i=0; i<Pxyz*Vt*sizeof(Gauge)/sizeof(fptype); ++i)
      if((i/(sizeof(Gauge)/sizeof(fptype)/8)) & 1) g_norm += ((fptype*)gfl[j])[i] * ((fptype*)gfl[j])[i];
    printf("gfl_norm[%d] = %f\n", j, g_norm);
  }

  for(int j=0; j<2; ++j)
  { g_norm = 0.0;
    for(int i=0; i<Pxyz*Vt*sizeof(Gauge)/sizeof(fptype); ++i)      if((i/(sizeof(Gauge)/sizeof(fptype)/8)+1) & 1)
      g_norm += ((fptype*)gfl[j])[i] * ((fptype*)gfl[j])[i];
    printf("gflback_norm[%d] = %f\n", j, g_norm);
  }
*/

  /* initialization process */
 start:

  int num_offsets_now = num_offsets;

//  printf("calling calc_source_norm...\n"); fflush(stdout);
#if TIME_CG
  t_misc = hrtimer_get();
#endif
  source_norm = calc_source_norm((fptype*)t_src);
#if TIME_CG
  //cout << "calc_source_norm time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif

  //printf("source_norm = %f\n", source_norm);

  for(int j=0; j<num_offsets; ++j) for(int i=0; i<num_floats_in_ks_array; ++i)
    ((fptype*)t_dest[j])[i] = 0.0;

/*
  for(int i=0; i<num_floats_in_ks_array; ++i)
    printf("t_dest0 = %f\n", ((fptype*)t_dest[1])[i]);
*/
/*
  FORSOMEPARITY_OMP(i,s,l_parity,private(j) reduction(+:source_norm) ){
    source_norm += (double) magsq_su3vec( src+i );
    su3vec_copy( src+i, resid+i);
    su3vec_copy(resid+i, cg_p+i);
    for(j=0;j<num_offsets;j++) {
	clearvec(psim[j]+i);
	su3vec_copy(resid+i, pm[j]+i);
      }
  } END_LOOP_OMP;
*/
  //copyKSn((fptype **)pm, (fptype *)t_src, num_offsets);
#if TIME_CG
  t_misc = hrtimer_get();
#endif
  for(j=0; j<num_offsets; ++j)
	copyKS((fptype *)pm[j], (fptype *)t_src);
  copyKS((fptype *)resid, (fptype *)t_src);
#if TIME_CG
  //cout << "copyKS time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif

/*
  for(int i=0; i<num_floats_in_ks_array; ++i) 
    printf("pm[%d] = %f\n", i, ((fptype*)pm[j_low])[i]);

  int padinc=0;
#if QPHIX_PrecisionInt==1
  for(int vt=0; vt<2; ++vt)
#else
  for(int vt=0; vt<1;  ++vt)
#endif 
    for(int tt=0; tt<Vt; ++tt) 
  for(int vz=0; vz<2; ++vz)  
    for(int zz=0; zz<Vz; ++zz)
  {
  for(int vy=0; vy<2; ++vy) 
    for(int yy=0; yy<Vy; ++yy)
  for(int vx=0; vx<2; ++vx)
    for(int xx=0; xx<Vxh; ++xx)
      for(int dir=0; dir<4; ++dir) for(int r=0; r<3; ++r) for(int c=0; c<3; ++c) 
    {
	padinc = tt*Vz+zz;
	printf("glle[%d][%d][%d][%d] = %f + I * %f\n", xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), dir, r, c, gll[0][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2+1][r][c][0][vx+2*vy+4*vz+8*vt], gll[0][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2+1][r][c][1][vx+2*vy+4*vz+8*vt]);
//	printf("gfle[%d][%d][%d][%d] = %f + I * %f\n", xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), dir, r, c, gfl[0][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2+1][r][c][0][vx+2*vy+4*vz+8*vt], gfl[0][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2+1][r][c][1][vx+2*vy+4*vz+8*vt]);
    }
//    if(vz==0&&vt==0) padinc++;
  }

  padinc=0;
#if QPHIX_PrecisionInt==1
  for(int vt=0; vt<2; ++vt)
#else
  for(int vt=0; vt<1; ++vt)
#endif
    for(int tt=0; tt<Vt; ++tt)
  for(int vz=0; vz<2; ++vz)
    for(int zz=0; zz<Vz; ++zz)
  {
  for(int vy=0; vy<2; ++vy)
    for(int yy=0; yy<Vy; ++yy)
  for(int vx=0; vx<2; ++vx)
    for(int xx=0; xx<Vxh; ++xx)
      for(int dir=0; dir<4; ++dir) for(int r=0; r<3; ++r) for(int c=0; c<3; ++c)
    {
	padinc = tt*Vz+zz;
	printf("glleback[%d][%d][%d][%d] = %f + I * %f\n", xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), dir, r, c, gll[0][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2][r][c][0][vx+2*vy+4*vz+8*vt], gll[0][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2][r][c][1][vx+2*vy+4*vz+8*vt]);
//	printf("gfleback[%d][%d][%d][%d] = %f + I * %f\n", xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), dir, r, c, gfl[0][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2][r][c][0][vx+2*vy+4*vz+8*vt], gfl[0][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2][r][c][1][vx+2*vy+4*vz+8*vt]);
    }
  //  if(vz==0&&vt==0) padinc++;
  }

  padinc=0;
#if QPHIX_PrecisionInt==1
  for(int vt=0; vt<2; ++vt)
#else
  for(int vt=0; vt<1; ++vt)
#endif
    for(int tt=0; tt<Vt; ++tt)
  for(int vz=0; vz<2; ++vz)
    for(int zz=0; zz<Vz; ++zz)
  {
  for(int vy=0; vy<2; ++vy)
    for(int yy=0; yy<Vy; ++yy)
  for(int vx=0; vx<2; ++vx)
    for(int xx=0; xx<Vxh; ++xx)
      for(int dir=0; dir<4; ++dir) for(int r=0; r<3; ++r) for(int c=0; c<3; ++c)
    {
	padinc = tt*Vz+zz;
	printf("gllo[%d][%d][%d][%d] = %f + I * %f\n", xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), dir, r, c, gll[1][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2+1][r][c][0][vx+2*vy+4*vz+8*vt], gll[1][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2+1][r][c][1][vx+2*vy+4*vz+8*vt]);
//	printf("gflo[%d][%d][%d][%d] = %f + I * %f\n", xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), dir, r, c, gfl[1][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2+1][r][c][0][vx+2*vy+4*vz+8*vt], gfl[1][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2+1][r][c][1][vx+2*vy+4*vz+8*vt]);
    }
  //if(vz==0&&vt==0) padinc++;
  }

  padinc=0;
#if QPHIX_PrecisionInt==1
  for(int vt=0; vt<2; ++vt)
#else
  for(int vt=0; vt<1; ++vt)
#endif
    for(int tt=0; tt<Vt; ++tt)
  for(int vz=0; vz<2; ++vz)
    for(int zz=0; zz<Vz; ++zz)
  {
  for(int vy=0; vy<2; ++vy)
    for(int yy=0; yy<Vy; ++yy)
  for(int vx=0; vx<2; ++vx)
    for(int xx=0; xx<Vxh; ++xx)
      for(int dir=0; dir<4; ++dir) for(int r=0; r<3; ++r) for(int c=0; c<3; ++c)
    {
	padinc = tt*Vz+zz;
	printf("glloback[%d][%d][%d][%d] = %f + I * %f\n", xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), dir, r, c, gll[1][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2][r][c][0][vx+2*vy+4*vz+8*vt], gll[1][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2][r][c][1][vx+2*vy+4*vz+8*vt]);
//    if(vz==0&&vt==0) padinc++;
//	printf("gfloback[%d][%d][%d][%d] = %f + I * %f\n", xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), dir, r, c, gfl[1][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2][r][c][0][vx+2*vy+4*vz+8*vt], gfl[1][padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][dir*2][r][c][1][vx+2*vy+4*vz+8*vt]);
    }
  }
*/
/*
  for(int i=0; i<num_floats_in_ks_array; ++i)
    printf("resid0[%d] = %f\n", i, ((fptype *)resid)[i]);
*/
  g_doublesum( &source_norm );
  rsq = source_norm;
  
  /* Special case -- trivial solution */
  if(source_norm == 0.){
    for(j = 0; j < num_offsets; j++){
      res_arg[j]->final_rsq     = 0.;
      res_arg[j]->size_r        = 0.;
      res_arg[j]->size_relr     = 0.;
      res_arg[j]->final_iter   = iteration;
      res_arg[j]->final_restart = nrestart;
    }
   
    //printf("multicg_qphix: source_norm = 0.0 parity = %d\n", parity);

    gettimeofday( &tv, NULL );
    info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec - info->final_sec;
    info->final_flop = (1205. + 15.*num_offsets) * total_iters * qphix_even_sites_on_node;
    info->status = QPHIX_SUCCESS; 
    /* Free stuff */
#ifdef TIME_CG
  t_misc = hrtimer_get();
#endif
    for(j=0;j<num_offsets;j++) freeKS(pm[j]); 
    freeKS(ttt);
    freeKS(ttt_temp);
    freeKS(resid);
    free(zeta_i);
    free(zeta_ip1);
    free(zeta_im1);
    free(beta_i);
    free(beta_im1);
    free(alpha);
    free(shifts);
    free(finished);
#ifdef TIME_CG
  //cout << "Free time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif
/*
  delete gBar;
  delete phaser;
  gBar=0x00;
  phaser=0x00;
*/
#if TIME_CG
    extern int myRank;
    if(myRank == 0)
    {
	double cgt = (hrtimer_get() - t)/1.0e09;
        std::cout << myname << ": CG_time(s) = " << cgt << " , Gflops = " << info->final_flop/1.e09/cgt << "\n";
    }
#if BW_CG       
    print_cg_timings(cgPerfData, num_offsets);
#endif
#endif
    return (iteration);
  }
    
  iteration++ ;  /* iteration counts number of multiplications
		    by M_adjoint*M */
#if BW_CG
  cgPerfData.iters++;
#endif

  total_iters++;
  double rsqstop = rsqmin * source_norm;
  
  for(j=0;j<num_offsets;j++){
    zeta_im1[j] = zeta_i[j] = 1.0;
    beta_im1[j] = -1.0;
    alpha[j] = 0.0;
  }
  
  do{
        oldrsq = rsq;
    /* sum of neighbors */
#if TIME_CG
        cgPerfData.dslash_time -= hrtimer_get();
#endif

/*
	for(int i=0; i<Pxyz*Vt; ++i)
	for(int dir=0; dir<8; ++dir) for(int c=0; c<3; ++c) for(int r=0; r<3; ++r) for(int v=0; v<VECLEN; ++v)
	  printf("gll[%d][%d][%d][%d][%d] = %f + I * %f\n", i, dir, c, r, v, gll[otherparity][i][dir][c][r][0][v], gll[otherparity][i][dir][c][r][1][v]);

        for(int i=0; i<Pxyz*Vt; ++i)
        for(int dir=0; dir<8; ++dir) for(int c=0; c<3; ++c) for(int r=0; r<3; ++r) for(int v=0; v<VECLEN; ++v)
          printf("gfl[%d][%d][%d][%d][%d] = %f + I * %f\n", i, dir, c, r, v, gfl[otherparity][i][dir][c][r][0][v], gfl[otherparity][i][dir][c][r][1][v]);

	printf("calling qphix_ks_dslash\n");
	fflush(stdout);
*/
        qphix_ks_dslash((void *) pm[j_low]
                      , (void *) gll[otherparity]
                      , (void *) gfl[otherparity]
                      , (void *) ttt_temp
                      , otherparity);

/*
  //if(parity0==QPHIX_ODD) 
  {
  int padinc=0;
#if QPHIX_PrecisionInt==1
  for(int vt=0; vt<2; ++vt)
#else
  for(int vt=0; vt<1; ++vt)
#endif
    for(int tt=0; tt<Vt; ++tt)
  for(int vz=0; vz<2; ++vz)
    for(int zz=0; zz<Vz; ++zz)
    {
	padinc = Vz*tt+zz;
    for(int vy=0; vy<2; ++vy)
    for(int yy=0; yy<Vy; ++yy)
    for(int vx=0; vx<2; ++vx)
    for(int xx=0; xx<Vxh; ++xx)
	for(int c=0; c<3; ++c)
	  printf("ttt_temp[%d][%d] = %f + I * %f\n",  xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), c, ttt_temp[padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][c][0][vx+2*vy+4*vz+8*vt], ttt_temp[padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][c][1][vx+2*vy+4*vz+8*vt]);
    //for(int c=0; c<3; ++c)
//	printf("ttt_temp[%d][%d] = %f + I * %f\n",  padinc+xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), c, ttt_temp[padinc+xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt))))))][c][0][vx+2*vy+4*vz+8*vt], ttt_temp[padinc+xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt))))))][c][1][vx+2*vy+4*vz+8*vt]);
    }
  }
*/
        qphix_ks_dslash((void *) ttt_temp
                      , (void *) gll[parity]
                      , (void *) gfl[parity]
                      , (void *) ttt
                      , parity);

/*
	for(int i=0; i<num_floats_in_ks_array; ++i)
	  printf("ttt[%d] = %f\n", i, ((fptype*)ttt)[i]);
*/
/*
  if(parity0==QPHIX_ODD) {
  int padinc=0;
  for(int vt=0; vt<2; ++vt)
    for(int tt=0; tt<Vt; ++tt)
  for(int vz=0; vz<2; ++vz)
    for(int zz=0; zz<Vz; ++zz)
    {
    for(int vy=0; vy<2; ++vy)
    for(int yy=0; yy<Vy; ++yy)
    for(int vx=0; vx<2; ++vx)
    for(int xx=0; xx<Vxh; ++xx)
        for(int c=0; c<3; ++c)
          printf("ttt[%d][%d] = %f + I * %f\n",  xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), c, ttt[padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][c][0][vx+2*vy+4*vz+8*vt], ttt[padinc+xx+Vxh*(yy+Vy*(zz+Vz*tt))][c][1][vx+2*vy+4*vz+8*vt]);
    if(vz==0&&vt==0)
    {
    //for(int c=0; c<3; ++c)
      //  printf("ttt[%d][%d] = %f + I * %f\n",  padinc+xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt)))))), c, ttt[padinc+xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt))))))][c][0][vx+2*vy+4*vz+8*vt], ttt[padinc+xx+Vxh*(vx+2*(yy+Vy*(vy+2*(zz+Vz*(vz+2*(tt+Vt*vt))))))][c][1][vx+2*vy+4*vz+8*vt]);
    padinc++;
    }
    }
  }
*/
#if TIME_CG
        cgPerfData.dslash_time += hrtimer_get();
#endif

/*    
	for(int i=0; i<num_floats_in_ks_array; ++i)
	  printf("ttt[%d] = %f\n", i, ttt[i]);
*/  
    /* finish computation of (-1)*M_adjoint*m*p and (-1)*p*M_adjoint*M*p */
    /* ttt  <- ttt - shift0*cg_p	*/
    /* pkp  <- cg_p . ttt */
/*
    pkp = 0.0;
    FORSOMEPARITY_OMP(i,s,l_parity,reduction(+:pkp) ){
      scalar_mult_add_su3_vector( ttt+i, cg_p+i, shift0, ttt+i );
      pkp += (double)su3_rdot( cg_p+i, ttt+i );
    } END_LOOP_OMP;
*/

#if TIME_CG
        cgPerfData.calc_pkp_time -= hrtimer_get();
#endif
        /* pkp  <- cg_p.(ttt - msq*cg_p) */
        pkp = calc_pkp((fptype*)ttt, (fptype*)pm[j_low], -1.0*offset_low);

/*
	printf("pkp = %f\n", pkp);
	fflush(stdout);
*/
#if TIME_CG
        cgPerfData.calc_pkp_time += hrtimer_get();
#endif
    g_doublesum( &pkp );
    iteration++;
    total_iters++;
#if BW_CG
    cgPerfData.iters++;
#endif

    beta_i[j_low] = -rsq / pkp;
    zeta_ip1[j_low] = 1.0; // this doesn't change (j_low vector is ordinary one mass CG)
    for(j=0;j<num_offsets_now;j++) if(j!=j_low){
	zeta_ip1[j] = zeta_i[j] * zeta_im1[j] * beta_im1[j_low];
	fptype c1 = beta_i[j_low] * alpha[j_low] * (zeta_im1[j]-zeta_i[j]);
	fptype c2 = zeta_im1[j] * beta_im1[j_low] * (1.0+shifts[j]*beta_i[j_low]);
	/*THISBLOWSUP
	  zeta_ip1[j] /= c1 + c2;
	  beta_i[j] = beta_i[j_low] * zeta_ip1[j] / zeta_i[j];
	*/
	/*TRYTHIS*/
	if( c1+c2 != 0.0 )
	  zeta_ip1[j] /= c1 + c2; 
	else {
	  zeta_ip1[j] = 0.0;
	  finished[j] = 1;
	}
	if( zeta_i[j] != 0.0){
	  beta_i[j] = beta_i[j_low] * zeta_ip1[j] / zeta_i[j];
	} else  {
	  zeta_ip1[j] = 0.0;
	  beta_i[j] = 0.0;
	  finished[j] = 1;
	}
	if(finished[j]==1&&j==num_offsets_now-1 && j>j_low) num_offsets_now--;
	//printf("beta_i[%d] = %f\n zeta_i[%d] = %f\n zeta_ip1[%d] = %f\n", j, beta_i[j], j, zeta_i[j], j, zeta_ip1[j]);
      }
    
    /* dest <- dest + beta*cg_p ( cg_p = pm[j_low], dest = psim[j_low] ) */
/*
    rsq = 0.0;
    FORSOMEPARITY_OMP(i,s,l_parity,private(j) reduction(+:rsq) ){
      for(j=0;j<num_offsets_now;j++) {
	scalar_mult_add_su3_vector( psim[j]+i, pm[j]+i, (Real)beta_i[j], psim[j]+i);
      }

      scalar_mult_add_su3_vector( resid+i, ttt+i, (Real)beta_i[j_low], resid+i);
      rsq += (double)magsq_su3vec( resid+i );
    } END_LOOP_OMP;
*/
#if TIME_CG
        cgPerfData.calc_rsq2_time -= hrtimer_get();
#endif

/*
  for(int i=0; i<num_floats_in_ks_array; ++i)
    printf("resid0[%d] = %f\n", i, ((fptype *)resid)[i]);
*/
	/* resid <- resid + a*ttt */
        rsq = calc_rsq2( (fptype*)resid
                           , beta_i[j_low]
                           , (fptype*)ttt
                       );

/*
	printf("rsq = %f\n", rsq);
	fflush(stdout);
*/
/*
	for(int i=0; i<num_floats_in_ks_array; ++i)
	  printf("resid[%d] = %f\n", i, ((fptype*)resid)[i]);

	exit(0);
*/
#if TIME_CG
        cgPerfData.calc_rsq2_time += hrtimer_get();
#endif

    g_doublesum(&rsq);
    
    if( rsq <= rsqstop ){

      /* calculate of t_dest */
      /* dest <- dest + a*cg_p */
#if TIME_CG
  t_misc = hrtimer_get();
#endif
      update_dest( t_dest
	, num_offsets_now
	, (const KS **)pm
	, beta_i
      );	
#if TIME_CG
  //cout << "update_dest time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif

      /* if parity==QPHIX_EVENODD, set up to do odd sites and go back */
      if(parity0 == QPHIX_EVENODD) {
	otherparity = parity0 & 1;
	parity = 1 - otherparity;
	parity0 = QPHIX_EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	t_src = (KS*)in->odd;
	for(j=0; j<num_offsets; ++j)
	    t_dest[j] = (KS*)out[j]->odd;
	goto start;
      }
      for(j = 0; j < num_offsets; j++){
	res_arg[j]->final_rsq     = (fptype)rsq/source_norm;
	res_arg[j]->size_r        = res_arg[j]->final_rsq;
	res_arg[j]->size_relr     = res_arg[j]->final_rel;
	res_arg[j]->final_iter   = iteration;
	res_arg[j]->final_restart = nrestart;
      }
     
      gettimeofday( &tv, NULL );
      info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec - info->final_sec;
      info->final_flop = (1205. + 15.*num_offsets) * total_iters * qphix_even_sites_on_node;
      info->status = QPHIX_SUCCESS; 
      /* Free stuff */
#ifdef TIME_CG
  t_misc = hrtimer_get();
#endif
      for(j=0;j<num_offsets;j++) freeKS(pm[j]); 
      freeKS(ttt);
      freeKS(ttt_temp);
      freeKS(resid);
      free(zeta_i);
      free(zeta_ip1);
      free(zeta_im1);
      free(beta_i);
      free(beta_im1);
      free(alpha);
      free(shifts);
      free(finished);
#ifdef TIME_CG
  //cout << "Free time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif
/*     
  delete gBar;
  delete phaser;
  gBar=0x00;
  phaser=0x00;
*/
#if TIME_CG
    extern int myRank;
    if(myRank == 0)
    {
	double cgt = (hrtimer_get() - t)/1.0e09;
        std::cout << myname << ": CG_time(s) = " << cgt << " , Gflops = " << info->final_flop/1.e09/cgt << "\n";
    }
/*
    if(myRank == 0) 
	std::cout << myname << ": Total CG Time (s) : " << (hrtimer_get() - t)/1.0e09 << "\n";
*/
#if BW_CG       
    print_cg_timings(cgPerfData, num_offsets);
#endif
#endif
/*
    for(int j=0; j<num_offsets; ++j) for(int i=0; i<num_floats_in_ks_array; ++i)
	printf("out[%d][%d][%d] = %f\n", parity, j, i, ((fptype*)t_dest[j])[i]);
*/
/*
    for(int jj=0; jj<num_offsets; ++jj) {
	double check = 0.0;
 	for(int i=0; i<num_floats_in_ks_array; ++i)
	  check += ((fptype*)t_dest[jj])[i] * ((fptype*)t_dest[jj])[i];
    	if(check==0.0) printf("t_dest[%d] = 0.0 for parity %d\n", jj, inv_arg->parity);
    }
*/
      return (iteration);
    }
    
    alpha[j_low] = rsq / oldrsq;

    //printf("alpha[%d] = %f\n", j_low, alpha[j_low]);
    
    for(j=0;j<num_offsets_now;j++) if(j!=j_low){
	/*THISBLOWSUP
	  alpha[j] = alpha[j_low] * zeta_ip1[j] * beta_i[j] /
	  (zeta_i[j] * beta_i[j_low]);
	*/
	/*TRYTHIS*/
	if( zeta_i[j] * beta_i[j_low] != 0.0)
	  alpha[j] = alpha[j_low] * zeta_ip1[j] * beta_i[j] /
	    (zeta_i[j] * beta_i[j_low]);
	else {
	  alpha[j] = 0.0;
	  finished[j] = 1;
	}
	//printf("alpha[%d] = %f\n", j, alpha[j]);
    }
    
    /* cg_p  <- resid + alpha*cg_p */
/*
    FORSOMEPARITY_OMP(i,s,l_parity,private(j) ){
      for(j=0;j<num_offsets_now;j++) {
	scalar_mult_su3_vector( resid+i, (Real)zeta_ip1[j], ttt+i);
	scalar_mult_add_su3_vector( ttt+i, pm[j]+i, (Real)alpha[j], pm[j]+i);
      }
      su3vec_copy(pm[j_low]+i,cg_p+i);
    } END_LOOP_OMP;
*/
#if TIME_CG
        cgPerfData.axpy_time -= hrtimer_get();
#endif
/*
	printf("calling update_pm\n");
	fflush(stdout);
*/
        update_pm( (const fptype *)resid
                 , t_dest /*, (fptype*)ttt*/
                 , pm
		 , beta_i
                 , zeta_ip1
                 , alpha
                 , num_offsets_now
                 );

#if TIME_CG
        cgPerfData.axpy_time += hrtimer_get();
#endif
    
    /* scroll the scalars */
    for(j=0;j<num_offsets_now;j++){
      beta_im1[j] = beta_i[j];
      zeta_im1[j] = zeta_i[j];
      zeta_i[j] = zeta_ip1[j];
    }
    
  } while( iteration < max_cg );
  
  for(j = 0; j < num_offsets; j++){
    res_arg[j]->final_rsq     = (fptype)rsq/source_norm;
    res_arg[j]->size_r        = res_arg[j]->final_rsq;
    res_arg[j]->size_relr     = res_arg[j]->final_rel;
    res_arg[j]->final_iter   = iteration;
    res_arg[j]->final_restart = nrestart;
    //->converged     = 0;
  }

#ifdef TIME_CG  
  extern int myRank;
  if(myRank == 0) 
	printf(
	       "%s: CG NOT CONVERGED after %d iterations, res. = %e wanted %e\n",
	       myname, iteration, res_arg[0]->final_rsq, rsqmin);
  fflush(stdout);
#endif
  
  /* if parity==QPHIX_EVENODD, set up to do odd sites and go back */
  if(parity0 == QPHIX_EVENODD) {
    otherparity = parity0 & 1;
    parity = 1 - otherparity;
    parity0 = QPHIX_EVEN;	/* so we won't loop endlessly */
    t_src = (KS*)in->odd;
    for(j=0; j<num_offsets; ++j)
	t_dest[j] = (KS*)out[j]->odd;
    iteration = 0;
    goto start;
  }
 
  gettimeofday( &tv, NULL );
  info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec - info->final_sec;
  info->final_flop = (1205. + 15.*num_offsets) * total_iters * qphix_even_sites_on_node;
  info->status = QPHIX_SUCCESS; 
  /* Free stuff */
#ifdef TIME_CG
  t_misc = hrtimer_get();
#endif
  for(j=0;j<num_offsets;j++){ freeKS(pm[j]); }
  freeKS(ttt);
  freeKS(ttt_temp);
  freeKS(resid);
  free(zeta_i);
  free(zeta_ip1);
  free(zeta_im1);
  free(beta_i);
  free(beta_im1);
  free(alpha);
  free(shifts);
  free(finished);
#ifdef TIME_CG
  //cout << "Free time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif
/*
  delete gBar;
  delete phaser;
  gBar=0x00;
  phaser=0x00;
*/
#if TIME_CG
  extern int myRank;
    if(myRank == 0)
    {
	double cgt = (hrtimer_get() - t)/1.0e09;
        std::cout << myname << ": CG_time(s) = " << cgt << ", Gflops = " << info->final_flop/1.e09/cgt << "\n";
    }
/*
  if(myRank == 0)
  {
        double cg_time = (double)(hrtimer_get() - t)/1.0e09;
        std::cout << myname << ": Total CG Time (s) : " << cg_time << " , Gflops : " << info->final_flop/cg_time << "\n";
  }
*/
#if BW_CG
    print_cg_timings(cgPerfData, num_offsets);
#endif
#endif  
/*
    for(int j=0; j<num_offsets; ++j) for(int i=0; i<num_floats_in_ks_array; ++i)
        printf("out[%d][%d][%d] = %f\n", parity, j, i, ((fptype*)t_dest[j])[i]);
*/
/*
    for(int j=0; j<num_offsets; ++j) {
        double check = 0.0;
        for(int i=0; i<num_floats_in_ks_array; ++i)
          check += ((fptype*)t_dest[j])[i] * ((fptype*)t_dest[j])[i];
        if(check==0.0) printf("t_dest[%d] = 0.0 for parity %d\n", j, inv_arg->parity);
  }
*/
  return(iteration);
}

