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
#define CG_DEBUG 0

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

extern int qphix_sites_on_node;

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

/*
    shifts = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_i = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_im1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_ip1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    beta_i = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    beta_im1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    alpha = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
*/
    shifts = (fptype * )malloc(num_offsets*sizeof(fptype));
    zeta_i = (fptype * )malloc(num_offsets*sizeof(fptype));
    zeta_im1 = (fptype * )malloc(num_offsets*sizeof(fptype));
    zeta_ip1 = (fptype * )malloc(num_offsets*sizeof(fptype));
    beta_i = (fptype * )malloc(num_offsets*sizeof(fptype));
    beta_im1 = (fptype * )malloc(num_offsets*sizeof(fptype));
    alpha = (fptype * )malloc(num_offsets*sizeof(fptype));

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

  /* initialization process */
 start:

  int num_offsets_now = num_offsets;

#if TIME_CG
  t_misc = hrtimer_get();
#endif
  source_norm = calc_source_norm((fptype*)t_src);
#if TIME_CG
  //cout << "calc_source_norm time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif

  for(int j=0; j<num_offsets; ++j) for(int i=0; i<num_floats_in_ks_array; ++i)
    ((fptype*)t_dest[j])[i] = 0.0;

#if TIME_CG
  t_misc = hrtimer_get();
#endif
  for(j=0; j<num_offsets; ++j)
	copyKS((fptype *)pm[j], (fptype *)t_src);
  copyKS((fptype *)resid, (fptype *)t_src);
#if TIME_CG
  //cout << "copyKS time : " << (double)(hrtimer_get() - t_misc) << "\n";
#endif

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
   
    gettimeofday( &tv, NULL );
    info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec - info->final_sec;
    info->final_flop = (1205. + 15.*num_offsets) * total_iters * qphix_sites_on_node;
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
  
/*  
  iteration++ ;
#if BW_CG
  cgPerfData.iters++;
#endif

  total_iters++;
*/
  double rsqstop = rsqmin * source_norm;
#ifdef CG_DEBUG
  if(myRank == 0)printf("%s: source_norm = %e\n", myname, (double)source_norm);
  fflush(stdout);
#endif
  
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

        qphix_ks_dslash((void *) pm[j_low]
                      , (void *) gll[otherparity]
                      , (void *) gfl[otherparity]
                      , (void *) ttt_temp
                      , otherparity);

        qphix_ks_dslash((void *) ttt_temp
                      , (void *) gll[parity]
                      , (void *) gfl[parity]
                      , (void *) ttt
                      , parity);

#if TIME_CG
        cgPerfData.dslash_time += hrtimer_get();
#endif

    /* finish computation of (-1)*M_adjoint*m*p and (-1)*p*M_adjoint*M*p */
    /* ttt  <- ttt - shift0*cg_p	*/
    /* pkp  <- cg_p . ttt */
#if TIME_CG
        cgPerfData.calc_pkp_time -= hrtimer_get();
#endif
        /* pkp  <- cg_p.(ttt - msq*cg_p) */
        pkp = calc_pkp((fptype*)ttt, (fptype*)pm[j_low], -1.0*offset_low);

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
	//	printf("beta_i[%d] = %f\n zeta_i[%d] = %f\n zeta_ip1[%d] = %f\n", j, beta_i[j], j, zeta_i[j], j, zeta_ip1[j]);
	//	fflush(stdout);
      }
    
    /* dest <- dest + beta*cg_p ( cg_p = pm[j_low], dest = psim[j_low] ) */
#if TIME_CG
        cgPerfData.calc_rsq2_time -= hrtimer_get();
#endif

	/* resid <- resid + a*ttt */
        rsq = calc_rsq2( (fptype*)resid
                           , beta_i[j_low]
                           , (fptype*)ttt
                       );

#if TIME_CG
        cgPerfData.calc_rsq2_time += hrtimer_get();
#endif

    g_doublesum(&rsq);
    
#if CG_DEBUG
    if(myRank==0)printf("%s: iter %d rsq = %g offsets = %d\n", 
			myname, iteration, rsq, num_offsets_now);
    fflush(stdout);
#endif
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
      info->final_flop = (1205. + 15.*num_offsets) * total_iters * qphix_sites_on_node;
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
      return (iteration);
    }
    
    alpha[j_low] = rsq / oldrsq;

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
#if TIME_CG
        cgPerfData.axpy_time -= hrtimer_get();
#endif
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
  info->final_flop = (1205. + 15.*num_offsets) * total_iters * qphix_sites_on_node;
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

#if TIME_CG
  extern int myRank;
    if(myRank == 0)
    {
	double cgt = (hrtimer_get() - t)/1.0e09;
        std::cout << myname << ": CG_time(s) = " << cgt << ", Gflops = " << info->final_flop/1.e09/cgt << "\n";
    }
#if BW_CG
    print_cg_timings(cgPerfData, num_offsets);
#endif
#endif  

  return(iteration);
}

