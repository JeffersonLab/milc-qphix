#include <string.h>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <cassert>
#include <time.h>

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

extern int qphix_sites_on_node;

/******************************** BLAS Kernels ********************************/
static double
calc_source_norm (fptype *restrict spinor)
{
    fptype source_norm = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:source_norm)    
    for(int i = 0; i < num_floats_in_ks_array; ++i) {
        source_norm += spinor[i]*spinor[i];    
    }
    return double(source_norm);
}

/* Original MILC code 
FORSOMEPARITYDOMAIN(i,s,parity){
    scalar_mult_add_su3_vector(&ttt[i],&t_dest[i],-msq_x4,&ttt[i]);
    add_su3_vector( &t_src[i], &ttt[i], &resid[i] );
    // remember ttt contains -M_adjoint*M*src 
    cg_p[i] = resid[i];
    rsq += (double)magsq_su3vec( &resid[i] );
} 
*/
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

static double
calc_rsq0 ( fptype *restrict ttt
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
        rsq += cg_p[i] * cg_p[i];
    }
    return double(rsq);
}

/* Original MILC code
FORSOMEPARITYDOMAIN(i,s,parity){
    scalar_mult_add_su3_vector( &t_dest[i], &cg_p[i], a, &t_dest[i] );
    scalar_mult_add_su3_vector( &resid[i], &ttt[i], a, &resid[i]);
#ifdef FEWSUMS
    actual_rsq += (double)magsq_su3vec( &resid[i] );
#else
    rsq += (double)magsq_su3vec( &resid[i] );
#endif
}
*/
static double
calc_rsq2 (fptype *restrict t_dest
          , const fptype *restrict cg_p
          , fptype a
          , fptype *restrict resid
          , const fptype *restrict ttt
          )
{
    fptype rsq = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:rsq)    
    for(int i = 0; i < num_floats_in_ks_array; ++i) {
        scalar_mult_add_su3_vector2(t_dest[i], cg_p[i], a, t_dest[i]);   
        scalar_mult_add_su3_vector2(resid[i], ttt[i], a, resid[i]);   
        rsq += resid[i] * resid[i];
    }
    return (double) rsq;
}

/* Original MILC code 
FORSOMEPARITYDOMAIN(i,s,parity) {
    scalar_mult_add_su3_vector( &ttt[i], &cg_p[i], -msq_x4, &ttt[i] );
    pkp += (double)su3_rdot( &cg_p[i], &ttt[i] );
#ifdef FEWSUMS
    c_tr += (double)su3_rdot( &ttt[i], &resid[i] );
    c_tt += (double)su3_rdot( &ttt[i], &ttt[i] );
#endif
}
*/
static double
calc_pkp ( fptype *restrict ttt
         , fptype *restrict cg_p
         , fptype msq_x4)
{
    fptype pkp = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:pkp)    
    for(int i = 0; i < num_floats_in_ks_array; ++i) {
        scalar_mult_add_su3_vector2(ttt[i], cg_p[i], -msq_x4, ttt[i]); 
        pkp += cg_p[i] * ttt[i];
    }
    return double(pkp);
}

/* Original MILC code        
FORSOMEPARITY(i,s,parity){
    scalar_mult_add_su3_vector( &resid[i], &cg_p[i], b, &cg_p[i]);
}
*/
static void
axpy ( const fptype *restrict resid
     , const fptype * restrict cg_p
     , fptype b
     , fptype *restrict ret)
{
#pragma omp parallel for num_threads (nThreads)
      for(int i = 0; i < num_floats_in_ks_array; ++i) {
        scalar_mult_add_su3_vector2(resid[i], cg_p[i], b, ret[i]); 
    }
}

/******************************************************************************/

#undef QPHIX_fptype_asqtad_invert
#undef QPHIX_fptype_FermionLinksAsqtad
#undef QPHIX_fptype_Real
#undef QPHIX_fptype_ColorVector

#if QPHIX_PrecisionInt==1
#define QPHIX_fptype_asqtad_invert QPHIX_F3_asqtad_invert
#define QPHIX_fptype_FermionLinksAsqtad QPHIX_F3_FermionLinksAsqtad
#define QPHIX_fptype_Real QPHIX_F_Real
#define QPHIX_fptype_ColorVector QPHIX_F3_ColorVector
#elif QPHIX_PrecisionInt==2
#define QPHIX_fptype_asqtad_invert QPHIX_D3_asqtad_invert
#define QPHIX_fptype_FermionLinksAsqtad QPHIX_D3_FermionLinksAsqtad
#define QPHIX_fptype_Real QPHIX_D_Real
#define QPHIX_fptype_ColorVector QPHIX_D3_ColorVector
#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif

int 
/*
qphix_ks_congrad_parity ( void *t_src_arg
                        , void *t_dest_arg
                        , struct QuarkInvertControl *qic
                        , fptype mass
                        , void *gll_arg[2]
                        , void *gfl_arg[2]
                        )
*/
QPHIX_fptype_asqtad_invert( QPHIX_info_t *info,
                        QPHIX_fptype_FermionLinksAsqtad *asqtad,
			QPHIX_invert_arg_t *inv_arg,
			QPHIX_resid_arg_t *res_arg,
			QPHIX_fptype_Real mass,
                        QPHIX_fptype_ColorVector *out,
                        QPHIX_fptype_ColorVector *in )
{
    if(myRank == 0 ) printf("QPHIX_fptype_asqtad_invert\n");
    fflush(stdout);
    char myname[25] = "qphix_ks_d_congrad_fn";

    assert(in->parity==out->parity && in->parity==inv_arg->parity);
#if TIME_CG
    hrtimer_t t = hrtimer_get();
    CGPerfData cgPerfData;
#endif

    Gauge *gll[2];
    Gauge18 *gfl[2];
    KS *t_src, *t_dest;
    if(in->parity!=QPHIX_ODD)
    {
	t_src = (KS*)in->even;
	t_dest = (KS*)out->even;
    }
    else
    {
	t_src = (KS*)in->odd;
	t_dest = (KS*)out->odd;
    }
    gll[0] = (Gauge*)asqtad->lngeven;
    gll[1] = (Gauge*)asqtad->lngodd ;
    
    gfl[0] = (Gauge18*)asqtad->fateven;
    gfl[1] = (Gauge18*)asqtad->fatodd;
    
    int iteration;                      /* counter for iterations */
    fptype a,b;                         /* Sugar's a,b */
#ifdef FEWSUMS
    double actual_rsq;                  /* rsq from actual summation of resid */
    double c_tr,c_tt,tempsum[4];        /* Re<resid|ttt>, <ttt|ttt> */
#endif
    double rsq = 0, relrsq = 0;         /* resid**2, rel resid*2 */
    double oldrsq, pkp;                 /* last resid*2, pkp = cg_p.K.cg_p */
    fptype msq_x4;                      /* 4*mass*mass */
    double source_norm;                 /* squared magnitude of source vector */
    int nrestart;                       /* Restart counter */

    /* temporary storages */
    KS *ttt, *ttt_temp, *cg_p, *resid;
    
    /* Unpack structure */
    int niter          = inv_arg->max;      
    int max_restarts   = inv_arg->nrestart; 
    fptype rsqmin      = res_arg->resid;    
    fptype relrsqmin   = res_arg->relresid; 
    int max_cg = max_restarts * niter;  /* Maximum number of iterations */
    msq_x4 = 4.0 * mass * mass;

    /* Convert MILC parity to QPhiX's parity. Unlike MILC, QPhiX only has 0 or 1 
     * as parity. */
    int parity         = inv_arg->parity & 1;   
    int otherparity    = 1 - parity;    /* the other parity */
    int parity0 = inv_arg->parity;

    /* Allocate the temps ttt, cg_p, resid */
    ttt      = (KS*)allocKS();
    ttt_temp = (KS*)allocKS();
    cg_p     = (KS*)allocKS();
    resid    = (KS*)allocKS();
  
    struct timeval tv;
    gettimeofday( &tv, NULL );
    /* Initialize info */
    info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec;
    info->final_flop = 0.;

    int total_iters = 0; 

    iteration = 0;

    start: 
    /* Source norm */
    source_norm = calc_source_norm((fptype*)t_src);
    g_doublesum( &source_norm );

    /* Start CG iterations */
    nrestart      = 0;
    total_iters += iteration;
    iteration     = 0;
    res_arg->size_r    = 0;
    res_arg->size_relr = 0;

    /* trivial source */
    if(source_norm ==0.0)
    {
	for(int i=0; i<num_floats_in_ks_array; ++i) 
	  ((fptype*)t_dest)[i] = 0.0;
	gettimeofday( &tv, NULL );	
	info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec - info->final_sec;

	//printf("d_congrad_qphix: source_norm = 0.0 parity = %d\n", parity);
#if TIME_CG
        extern int myRank;
        if(myRank == 0) std::cout << "Total CG Time (s) : " << (hrtimer_get() - t)/1.0e09 << "\n";
#endif	
	res_arg->final_rel = 0;
	res_arg->final_rsq = 0;
	res_arg->final_iter = total_iters;
	res_arg->final_restart = nrestart;

	/* Free */
	freeKS(ttt);
	freeKS(ttt_temp);
	freeKS(cg_p);
	freeKS(resid);

	return iteration;
    }

    while(true) {

        /* Check for completion */
        if(( iteration % niter == 0 ) || 
           (( rsqmin    <= 0 || rsqmin    > res_arg->size_r ) &&
            ( relrsqmin <= 0 || relrsqmin > res_arg->size_relr )
           )) 
        {
            /* (re)initialization process */
            /* Compute true residual and relative residual */
            rsq = 0.0;

#if TIME_CG
            cgPerfData.dslash_time -= hrtimer_get();
#endif            
	    //printf("calling qphix_ks_dslash ...\n");
	    //fflush(stdout);
            qphix_ks_dslash(t_dest
                          , gll[otherparity]
                          , gfl[otherparity]
                          , ttt_temp
                          , otherparity);
                          
            qphix_ks_dslash(ttt_temp
                          , gll[parity]
                          , gfl[parity]
                          , ttt
                          , parity);
                          
#if TIME_CG
            cgPerfData.dslash_time += hrtimer_get();
#endif   

#if TIME_CG
            cgPerfData.calc_rsq_time -= hrtimer_get();
#endif             
	    //printf("calling calc_rsq\n");
	    //fflush(stdout);
            rsq = calc_rsq( (fptype*)ttt
                          , (fptype*)t_dest
                          , (fptype*)t_src
                          , (fptype*)resid
                          , (fptype*)cg_p
                          , msq_x4
                          );
	    //printf("rsq = %f\n");
	    //fflush(stdout);
#if TIME_CG
            cgPerfData.calc_rsq_time += hrtimer_get();
#endif  
           
#ifdef FEWSUMS
            actual_rsq = rsq; /* not yet summed over nodes */
#endif
            g_doublesum( &rsq );
            if(relrsqmin > 0) {
#if TIME_CG
                cgPerfData.relative_residue_time -= hrtimer_get();
#endif              
                relrsq = relative_residue(resid, t_dest, parity);
#if TIME_CG
                cgPerfData.relative_residue_time += hrtimer_get();
#endif                  
            }
            res_arg->final_rsq    = (fptype)rsq/source_norm;
            res_arg->final_rel = (fptype)relrsq;

            /* iteration counts number of multiplications by M_adjoint*M */
            iteration++ ;
#if BW_CG
            cgPerfData.iters++;
#endif             
            //total_iters++;

            /* Quit when true residual and true relative residual are within
               tolerance or when we exhaust iterations or restarts */
            if(iteration >= max_cg || 
               nrestart  >= max_restarts ||
               ((rsqmin <= 0 || rsqmin > res_arg->final_rsq) &&
                (relrsqmin <= 0 || relrsqmin > res_arg->final_rel))) {

/*
#if TIME_CG
                extern int myRank;
                if(myRank == 0)
                {
                        double cg_time = (hrtimer_get() - t)/1.0e09;
			info->final_flop = 1187.*(1+total_iters)*qphix_sites_on_node;
                        std::cout << myname << ": Total CG Time (s) : " << cg_time << " , Gflops : " << info->final_flop/cg_time << "\n";
                }
#if BW_CG
                print_cg_timings(cgPerfData);
#endif
#endif
*/
		if(parity0==QPHIX_EVENODD) {
		    otherparity = parity0 & 1;
		    parity = 1 - otherparity;
		    parity0 = QPHIX_EVEN;
		    t_src = (KS*)in->odd;
		    t_dest = (KS*)out->odd;
		    goto start;
		}
		else break;
            }
            if(iteration!=0)  nrestart++;
        } // end of if (check for completion)

#ifdef FEWSUMS
        oldrsq = actual_rsq;    /* not yet summed over nodes */
#else
        oldrsq = rsq;
#endif
        /* sum of neighbours */

#if TIME_CG
        cgPerfData.dslash_time -= hrtimer_get();
#endif          
            //printf("calling qphix_ks_dslash ...\n");
            //fflush(stdout);
        qphix_ks_dslash(cg_p
                      , gll[otherparity]
                      , gfl[otherparity]
                      , ttt_temp
                      , otherparity);
                      
        qphix_ks_dslash(ttt_temp
                      , gll[parity]
                      , gfl[parity]
                      , ttt
                      , parity);        
#if TIME_CG
        cgPerfData.dslash_time += hrtimer_get();
#endif 
        /* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
        /* ttt  <- ttt - msq_x4*cg_p    (msq = mass squared) */

#ifdef FEWSUMS
        c_tr = 0.0; 
        c_tt = 0.0;
#endif

#if TIME_CG
        cgPerfData.calc_pkp_time -= hrtimer_get();
#endif         
        /* pkp  <- cg_p.(ttt - msq*cg_p) */
	//printf("calling calc_pkp\n");
	//fflush(stdout);
        pkp = calc_pkp((fptype*)ttt, (fptype*)cg_p, msq_x4);
	//printf("pkp = %f\n", pkp);
	//fflush(stdout);
#if TIME_CG
        cgPerfData.calc_pkp_time += hrtimer_get();
#endif 

#ifdef FEWSUMS
        /* finally sum oldrsq over nodes, also other sums */
        tempsum[0] = pkp; 
        tempsum[1] = c_tr; 
        tempsum[2] = c_tt; 
        tempsum[3] = oldrsq;
        g_vecdoublesum( tempsum, 4 );
        pkp = tempsum[0]; 
        c_tr = tempsum[1]; 
        c_tt = tempsum[2]; 
        oldrsq = tempsum[3];
#else
        g_doublesum( &pkp );
#endif
        iteration++;

#if BW_CG
        cgPerfData.iters++;
#endif                     
    
        a = (fptype) (-rsq/pkp);

#ifdef FEWSUMS
        actual_rsq = 0.0;
#else
        rsq = 0.0;
#endif

#if TIME_CG
        cgPerfData.calc_rsq2_time -= hrtimer_get();
#endif       
        /* dest <- dest + a*cg_p */  
        /* resid <- resid + a*ttt */
        rsq = calc_rsq2( (fptype*)t_dest
                       , (fptype*)cg_p
                       , a
                       , (fptype*)resid
                       , (fptype*)ttt
                       );
#if TIME_CG
        cgPerfData.calc_rsq2_time += hrtimer_get();
#endif        

#ifdef FEWSUMS
        /* TEST - should equal actual_rsq */
        rsq = oldrsq + 2.0*a*c_tr + a*a*c_tt; 
        /**c_tt = actual_rsq;**/ /* TEMP for test */
        /**g_doublesum(&c_tt);**/ /* TEMP true value for rsq */
#else
        g_doublesum(&rsq);
#endif    

        if(relrsqmin > 0) {
#if TIME_CG
            cgPerfData.relative_residue_time -= hrtimer_get();
#endif     
            relrsq = relative_residue(resid, t_dest, parity);
#if TIME_CG
            cgPerfData.relative_residue_time += hrtimer_get();
#endif       
        }
    
        res_arg->size_r    = (fptype)rsq/source_norm;
        res_arg->size_relr = (fptype)relrsq;
	res_arg->final_rsq     = res_arg->size_r;
        b = (fptype)rsq/oldrsq;

#if TIME_CG
        cgPerfData.axpy_time -= hrtimer_get();
#endif    
        /* cg_p  <- resid + b*cg_p */    
        axpy((fptype*)resid, (fptype*)cg_p, b, (fptype*)cg_p);

#if TIME_CG
        cgPerfData.axpy_time += hrtimer_get();
#endif    
    } // end of while

    gettimeofday( &tv, NULL ); 
    info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec - info->final_sec; 

    total_iters += iteration;

    res_arg->final_iter   = total_iters;
    res_arg->final_restart = nrestart;

<<<<<<< HEAD
    info->final_flop = 1187. * total_iters * qphix_even_sites_on_node;

=======
>>>>>>> gauge_force
    if(myRank==0) if(nrestart == max_restarts || iteration == max_cg) {
        if(omp_get_thread_num() == 0) {
            std::cout << "qphix_congrad: CG not converged after " << iteration
                      << " iterations and " << nrestart << " restarts, \n"
                      << "rsq. = " << std::scientific << res_arg->final_rsq 
                      << " wanted " << rsqmin << " relrsq = " 
                      << res_arg->final_rel << " wanted " << relrsqmin
                      << std::endl;
        }
        
    }
  
    /*! \fixme What is done here? Do we need it? */
    //cleanup_dslash_temps();
    freeKS(ttt); 
    freeKS(ttt_temp);
    freeKS(cg_p); 
    freeKS(resid);
/*
  delete gBar;
  delete phaser;
  gBar=0x00;
  phaser=0x00;
*/
/*
#if TIME_CG
		extern int myRank;
    if(myRank == 0) std::cout << "Total CG Time (s) : " << (hrtimer_get() - t)/1.0e09 << "\n";
#endif    
*/
#if TIME_CG
                extern int myRank;
                if(myRank == 0)
                {
                        double cg_time = (hrtimer_get() - t)/1.0e09;
                        info->final_flop = 1187.*total_iters*qphix_sites_on_node;
                        std::cout << myname << ": Total_CG_time(s) = " << cg_time << " , Gflops = " << info->final_flop/cg_time/1.0e09 << "\n";
                }
#if BW_CG
                print_cg_timings(cgPerfData);
#endif
#endif

    return iteration;
}



/*********** Added for multi-mass cg *************/

#if 0
static double
calc_rsq2_mul ( fptype *restrict t_dest[MAXNDMASS]
              , int num_offsets
              , int j_low
              , const fptype **restrict cg_p
              , fptype *a
              , fptype *restrict resid
              , const fptype *restrict ttt
              , int parity
              )
{
    fptype rsq = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:rsq) collapse(2) 
    for(int j=0; j<num_offsets; ++j) for(int i = 0; i < num_floats_in_ks_array; ++i) {
        //int pad = (parity==ODD ? num_floats_in_ks_array : 0);
        int pad = 0;
	scalar_mult_add_su3_vector2(t_dest[j][i+pad], cg_p[j][i], a[j], t_dest[j][i+pad]);
        if(j==j_low) {
	  if(i<50)
		printf("resid = %f\nttt = %f\n", resid[i], ttt[i]);
          scalar_mult_add_su3_vector2(resid[i], ttt[i], a[j_low], resid[i]);
	  if(i<50)
		printf("resid = %f\nttt = %f\n", resid[i], ttt[i]);
          rsq += resid[i] * resid[i];
        }
    }
    return (double) rsq;
}

static void
update_pm ( fptype *restrict resid
          , fptype *restrict ttt
          , fptype **restrict pm
          , fptype *restrict zeta_ip1
          , fptype *restrict alpha
          , int num_offsets
          )
{
#pragma omp parallel for num_threads (nThreads) collapse(2)
    for(int j=0; j<num_offsets; ++j) for(int i=0; i<num_floats_in_ks_array; ++i)
    {
        scalar_mult_su3_vector2(resid[i], zeta_ip1[j], ttt[i]);
        scalar_mult_add_su3_vector2(ttt[i], pm[j][i], alpha[j], pm[j][i]);
    }
}

int
qphix_ks_multicg_offset ( void *t_src_arg
                        , void **t_dest_arg
                        , struct QuarkInvertControl *qic
                        , fptype *mass
                        , int num_offsets
                        , void *gll_arg[2]
                        , void *gfl_arg[2]
                        )
{
#if TIME_CG
    hrtimer_t t = hrtimer_get();
    CGPerfData cgPerfData;
#endif
    Gauge *gll[2];
    Gauge18 *gfl[2];
    KS *t_src       = (KS*)t_src_arg;
    KS *t_dest[MAXNDMASS];
    for(int i=0; i<num_offsets; ++i)
    t_dest[i]      = (KS*)t_dest_arg[i];
    gll[0] = (Gauge*) gll_arg[0];
    gll[1] = (Gauge*) gll_arg[1];

    gfl[0] = (Gauge18*) gfl_arg[0];
    gfl[1] = (Gauge18*) gfl_arg[1];

    int iteration;                      /* counter for iterations */
    fptype a,b;                         /* Sugar's a,b */
#ifdef FEWSUMS
    double actual_rsq;                  /* rsq from actual summation of resid */
    double c_tr,c_tt,tempsum[4];        /* Re<resid|ttt>, <ttt|ttt> */
#endif
    double rsq = 0, relrsq = 0;         /* resid**2, rel resid*2 */
    double oldrsq, pkp;                 /* last resid*2, pkp = cg_p.K.cg_p */
    fptype msq_x4;                      /* 4*mass*mass */
    double source_norm;                 /* squared magnitude of source vector */
    int nrestart;                       /* Restart counter */

    /* temporary storages */
    KS *ttt, *ttt_temp, /**cg_p,*/ *resid;
    int j, j_low;
    fptype * restrict shifts, shift0;
    fptype * restrict zeta_i, * restrict zeta_im1, * restrict zeta_ip1;
    fptype * restrict beta_i, * restrict beta_im1, * restrict alpha;
    int * restrict finished;
    KS ** restrict pm;

    /* Unpack structure */
    int niter          = inv_arg->max;
    int max_restarts   = inv_arg->nrestart;
    fptype rsqmin      = resid_arg->resid;
    fptype relrsqmin   = resid_arg->relresid;
    int max_cg = max_restarts * niter;  /* Maximum number of iterations */

    //printf("niter = %d max_restarts = %d\n", niter, max_restarts);
    //printf("rsqmin = %.8e relrsqmin = %.8e\n", rsqmin, relrsqmin);

    /* Convert MILC parity to QPhiX's parity. Unlike MILC, QPhiX only has 0 or 1 
     * as parity. */

    int parity0 = inv_arg->parity, parity, otherparity; //        = inv_arg->parity & 1; 
    /* if we want both parities, we will do even first. */
    switch(parity0){
    case(QPHIX_EVEN): parity=QPHIX_EVEN; otherparity=QPHIX_ODD; break;
    case(QPHIX_ODD):  parity=QPHIX_ODD; otherparity=QPHIX_EVEN; break;
    case(QPHIX_EVENODD):  parity=QPHIX_EVEN; otherparity=QPHIX_ODD; break;
    }

    /* Allocate the temps ttt, cg_p, resid */
    ttt      = (KS*)allocKS();
    ttt_temp = (KS*)allocKS();
    resid    = (KS*)allocKS();

    shifts = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_i = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_im1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    zeta_ip1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    beta_i = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    beta_im1 = (fptype * restrict )malloc(num_offsets*sizeof(fptype));
    alpha = (fptype * restrict )malloc(num_offsets*sizeof(fptype));

    finished = (int *)malloc(sizeof(int)*num_offsets);
    /*
 *     if(finished == NULL){
 *           printf("%s: No room for 'finished'\n",myname);
 *                 terminate(1);
 *                     }
 *                         */
    double offset_low = 1.0e+20;
    j_low = -1;
    for(j=0;j<num_offsets;j++){
      assert((shifts[j] = mass[j])>=0.0);
      shifts[j] *= shifts[j] * 4.0;
      if (mass[j] < offset_low){
        offset_low = mass[j] * mass[j] * 4.0;
        j_low = j;
      }
    }
    msq_x4 = 4.0 * mass[j_low] * mass[j_low];

    for(j=0;j<num_offsets;j++){
      if( j!=j_low )shifts[j] -= shifts[j_low];
      }
    shift0 = -shifts[j_low];

    pm = (KS **)malloc(num_offsets*sizeof(KS *));
    for(j=0; j<num_offsets; ++j)
        pm[j] = (KS *)allocKS();

  start:

    for(j = 0; j < num_offsets; j++)
    {
        finished[j] = 0;
        beta_im1[j] = -1.;
        zeta_im1[j] = 1.;
        zeta_i[j] = 1.;
        alpha[j] = 0.;
    }
    /* Source norm */
    source_norm = calc_source_norm((fptype*)t_src /* + (parity == ODD ? num_floats_in_ks_array: 0)*/);
    g_doublesum( &source_norm );

    /* Start CG iterations */
    nrestart      = 0;
    iteration     = 0;
    int num_offsets_now = num_offsets;
    {
      resid_arg->size_r    = 0;
      resid_arg->size_relr = 0;
    }

    while(true) {
        /* Check for completion */
/*
 *         double size_r_max = resid_arg[0].size_r, size_relr_max = resid_arg[0].size_relr;
 *                 for(j = 1; j < num_offsets; j++){
 *                           if(size_r_max<resid_arg[j].size_r) size_r_max=resid_arg[j].size_r;
 *                                     if(size_relr_max<resid_arg[j].size_relr) size_relr_max=resid_arg[j].size_relr;
 *                                             }
 *                                             */
        if(( iteration % niter == 0 ) ||
           (( rsqmin    <= 0 || rsqmin    > resid_arg->size_r /*size_r_max*/ && resid_arg->size_r > 0.0) &&
            ( relrsqmin <= 0 || relrsqmin > resid_arg->size_relr /*size_relr_max*/ && resid_arg->size_relr > 0.0 )
           ))
        {
            /* (re)initialization process */
            /* Compute true residual and relative residual */
            rsq = 0.0;

#if TIME_CG
            cgPerfData.dslash_time -= hrtimer_get();
#endif
/*
	    for(j=0; j<num_offsets_now; ++j) for(int i=0; i<50; ++i) for(int c=0; c<3; ++c)
		printf("ks[%d][%d][%d] = %f + I * %f\n", j, i, c, t_dest[j][i][c][0][0], t_dest[j][i][c][1][0]);
            printf("calling qphix_ks_dslash ...\n");
*/
            qphix_ks_dslash(t_dest[j_low]
                          , gll[otherparity]
                          , gfl[otherparity]
                          , ttt_temp
                          , otherparity);

            qphix_ks_dslash(ttt_temp
                          , gll[parity]
                          , gfl[parity]
                          , ttt
                          , parity);

#if TIME_CG
            cgPerfData.dslash_time += hrtimer_get();
#endif

#if TIME_CG
            cgPerfData.calc_rsq_time -= hrtimer_get();
#endif
            rsq = calc_rsq( (fptype*)ttt
                          , (fptype*)t_dest[j_low]/* + (parity==ODD ? num_floats_in_ks_array : 0)*/
                          , (fptype*)t_src /*+ (parity==ODD ? num_floats_in_ks_array : 0)*/
                          , (fptype*)resid
                          , (fptype*)pm[j_low] /*cg_p*/
                          , msq_x4
                          );

	    //printf("rsq = %f\n", rsq);

#if TIME_CG
            cgPerfData.calc_rsq_time += hrtimer_get();
#endif

#ifdef FEWSUMS
            actual_rsq = rsq; /* not yet summed over nodes */
#endif
            g_doublesum( &rsq );
            if(relrsqmin > 0) {
#if TIME_CG
                cgPerfData.relative_residue_time -= hrtimer_get();
#endif
                relrsq = relative_residue(resid, t_dest[j_low] /*+ (parity==ODD ? num_floats_in_ks_array : 0)*/, parity);
#if TIME_CG
                cgPerfData.relative_residue_time += hrtimer_get();
#endif
            }

            /* iteration counts number of multiplications by M_adjoint*M */
  //          iteration++ ;

                resid_arg->final_rsq     = (fptype)rsq/source_norm;
                resid_arg->size_r        = resid_arg->final_rsq;
                resid_arg->size_relr     = resid_arg->final_relrsq;
                resid_arg->final_iters   = iteration;
/*
#if TIME_CG
            cgPerfData.iters++;
#endif
*/
            if(iteration >= max_cg ||
               nrestart  >= max_restarts ||
               ((rsqmin <= 0 || rsqmin > resid_arg->final_rsq) &&
                (relrsqmin <= 0 || relrsqmin > resid_arg->final_relrsq))) {

            	/* iteration counts number of multiplications by M_adjoint*M */
            	iteration++ ;
#if TIME_CG
                extern int myRank;
                if(myRank == 0)
                {
                        double cg_time = (hrtimer_get() - t)/1.0e09;
                        std::cout << myname << ": Total CG Time (s) : " << cg_time << " , Gflops : " << info->final_flop/cg_time << "\n";
                }
#if BW_CG
		cgPerfData.iters++;
                print_cg_timings(cgPerfData);
#endif
#endif
                if(parity0==QPHIX_EVENODD && rsqmin > resid_arg->final_rsq && relrsqmin > resid_arg->final_relrsq) {
                  parity=QPHIX_ODD;
                  otherparity=QPHIX_EVEN;
                  parity0=QPHIX_ODD;

                  goto start;

                }
                else
                  break;
            }
            nrestart++;
        } // end of if (check for completion)
#ifdef FEWSUMS
        oldrsq = actual_rsq;    /* not yet summed over nodes */
#else
        oldrsq = rsq;
#endif
        /* sum of neighbours */

#if TIME_CG
        cgPerfData.dslash_time -= hrtimer_get();
#endif

/*
	for(j=0; j<num_offsets_now; ++j) for(int i=0; i<50; ++i) for(int c=0; c<3; ++c)
	  printf("pm[%d][%d][%d] = %f + I * %f\n", j, i, c, pm[j][i][c][0][0], pm[j][i][c][1][0]);
	printf("calling qphix_ks_dslash ...\n");
*/
        qphix_ks_dslash(pm[j_low] /*cg_p*/
                      , gll[otherparity]
                      , gfl[otherparity]
                      , ttt_temp
                      , otherparity);

        qphix_ks_dslash(ttt_temp
                      , gll[parity]
                      , gfl[parity]
                      , ttt
                      , parity);
#if TIME_CG
        cgPerfData.dslash_time += hrtimer_get();
#endif

	start2:

#ifdef FEWSUMS
        c_tr = 0.0;
        c_tt = 0.0;
#endif

#if TIME_CG
        cgPerfData.calc_pkp_time -= hrtimer_get();
#endif

        /* pkp  <- cg_p.(ttt - msq*cg_p) */
        pkp = calc_pkp((fptype*)ttt, (fptype*)pm[j_low] /*cg_p*/, msq_x4);

	//printf("pkp = %f\n", pkp);

#if TIME_CG
        cgPerfData.calc_pkp_time += hrtimer_get();
#endif

#ifdef FEWSUMS
        /* finally sum oldrsq over nodes, also other sums */
        tempsum[0] = pkp;
        tempsum[1] = c_tr;
        tempsum[2] = c_tt;
        tempsum[3] = oldrsq;
        g_vecdoublesum( tempsum, 4 );
        pkp = tempsum[0];
        c_tr = tempsum[1];
        c_tt = tempsum[2];
        oldrsq = tempsum[3];
#else
        g_doublesum( &pkp );
#endif
        iteration++;

#if BW_CG
        cgPerfData.iters++;
#endif

        a = (fptype) (-rsq/pkp);

        beta_i[j_low] = (fptype) (-rsq/pkp);
        zeta_ip1[j_low] = 1.0;

        for(j=0;j<num_offsets_now;j++) if(j!=j_low){
          zeta_ip1[j] = zeta_i[j] * zeta_im1[j] * beta_im1[j_low];
          fptype c1 = beta_i[j_low] * alpha[j_low] * (zeta_im1[j]-zeta_i[j]);
          fptype c2 = zeta_im1[j] * beta_im1[j_low] * (1.0+shifts[j]*beta_i[j_low]);

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
            if(j==num_offsets_now-1 && j>j_low) num_offsets_now--;
          }
        }

#ifdef FEWSUMS
        actual_rsq = 0.0;
#else
        rsq = 0.0;
#endif

#if TIME_CG
        cgPerfData.calc_rsq2_time -= hrtimer_get();
#endif
        /* dest <- dest + a*cg_p */
        /* resid <- resid + a*ttt */
/*
	for(j=0; j<num_offsets_now; ++j)
	  printf("beta_i[%d] = %f\n", j, beta_i[j]);
*/
        rsq = calc_rsq2_mul( (fptype**)t_dest
                           , num_offsets_now
                           , j_low
                           , (const fptype**)pm /*cg_p*/
                           , beta_i
                           , (fptype*)resid
                           , (fptype*)ttt
                           , parity
                       );

        //printf("rsq_2 = %f\n", rsq);
#if TIME_CG
        cgPerfData.calc_rsq2_time += hrtimer_get();
#endif

#ifdef FEWSUMS
        /* TEST - should equal actual_rsq */
        rsq = oldrsq + 2.0*a*c_tr + a*a*c_tt;
        /**c_tt = actual_rsq;**/ /* TEMP for test */
        /**g_doublesum(&c_tt);**/ /* TEMP true value for rsq */
#else
        g_doublesum(&rsq);
#endif

        if(relrsqmin > 0) {
#if TIME_CG
            cgPerfData.relative_residue_time -= hrtimer_get();
#endif
            relrsq = relative_residue(resid, t_dest[j_low], parity);
#if TIME_CG
            cgPerfData.relative_residue_time += hrtimer_get();
#endif
        }

        resid_arg->size_r    = (fptype)rsq/source_norm;
        resid_arg->size_relr = (fptype)relrsq;

        alpha[j_low] = rsq / oldrsq;
        for(j=0;j<num_offsets_now;j++) if(j!=j_low){
          if( zeta_i[j] * beta_i[j_low] != 0.0)
            alpha[j] = alpha[j_low] * zeta_ip1[j] * beta_i[j] /
                        (zeta_i[j] * beta_i[j_low]);
          else {
            alpha[j] = 0.0;
            finished[j] = 1;
          }
        }

#if TIME_CG
        cgPerfData.axpy_time -= hrtimer_get();
#endif
        /* cg_p  <- resid + b*cg_p */
        update_pm( (fptype*)resid
                 , (fptype*)ttt
                 , (fptype**)pm
                 , zeta_ip1
                 , alpha
                 , num_offsets_now
                 );

/*
        for(j=0; j<num_offsets_now; ++j) for(int i=0; i<50; ++i) for(int c=0; c<3; ++c)
          printf("pm[%d][%d][%d] = %f + I * %f\n", j, i, c, pm[j][i][c][0][0], pm[j][i][c][1][0]);
*/
#if TIME_CG
        cgPerfData.axpy_time += hrtimer_get();
#endif

        /* scroll the scalars */
        for(j=0;j<num_offsets_now;j++){
          beta_im1[j] = beta_i[j];
          zeta_im1[j] = zeta_i[j];
          zeta_i[j] = zeta_ip1[j];

	  //printf("beta_i[%d] = %f zeta_i[%d] = %f zeta_ip1[%d] = %f\n", j, beta_i[j], j, zeta_i[j], j, zeta_ip1[j]);
        }

    } // end of while
    if(myRank==0) if(nrestart == max_restarts || iteration == max_cg) {
        if(omp_get_thread_num() == 0) {
            std::cout << "qphix_congrad: CG not converged after " << iteration
                      << " iterations and " << nrestart << " restarts, \n"
                      << "rsq. = " << std::scientific << resid_arg->final_rsq
                      << " wanted " << rsqmin << " relrsq = "
                      << resid_arg->final_relrsq << " wanted " << relrsqmin
                      << std::endl;
        }

    }

    /*! \fixme What is done here? Do we need it? */
    freeKS(ttt);
    freeKS(ttt_temp);
    freeKS(resid);
    for(j=0; j<num_offsets; ++j) 
	freeKS(pm[j]);
    free(pm);

#if TIME_CG
                extern int myRank;
    if(myRank == 0) std::cout << "Total CG Time (s) : " << (hrtimer_get() - t)/1.0e09 << "\n";
#endif
    return iteration;
}


static double
calc_rsq_mul ( fptype *restrict ttt
         , const fptype *restrict t_dest
         , const fptype *restrict t_src
         , fptype *restrict resid
         , fptype *restrict cg_p
         , fptype msq_x4)
{
    fptype rsq = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:rsq)
    for(int i = 0; i < num_floats_in_ks_array; ++i) {
        scalar_mult_add_su3_vector2(ttt[i], t_dest[i], msq_x4, ttt+i);
        add_su3_vector2(t_src[i], -ttt[i], resid+i);
        cg_p[i] = resid[i];
        rsq += resid[i] * resid[i];
    }
    return double(rsq);
}

static double
calc_rsq2_mul ( fptype *restrict t_dest[MAXNDMASS]
              , fptype *a
              , const fptype **restrict cg_p
              , fptype *restrict resid
              , const fptype *restrict ttt
              , int num_offsets
              , int j_low
              )
{
    fptype rsq = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:rsq) collapse(2) 
    for(int j=0; j<num_offsets; ++j) for(int i = 0; i < num_floats_in_ks_array; ++i) {
        int pad = 0;
        scalar_mult_add_su3_vector2(t_dest[j][i+pad], cg_p[j][i], -1.0*a[j], t_dest[j]+i+pad);
        if(j==j_low) {
          scalar_mult_add_su3_vector2(resid[i], ttt[i], a[j_low], resid[i]);
          rsq += resid[i] * resid[i];
        }
    }
    return (double) rsq;
}

static void
update_pm ( fptype **restrict pm
          , fptype *restrict zeta_ip1
          , fptype *restrict resid
          , fptype *restrict alpha
          , fptype *restrict ttt
          , int num_offsets
          )
{
#pragma omp parallel for num_threads (nThreads) collapse(2)
    for(int j=0; j<num_offsets; ++j) for(int i=0; i<num_floats_in_ks_array; ++i)
    {
        scalar_mult_su3_vector2(resid[i], zeta_ip1[j], ttt[i]);
        scalar_mult_add_su3_vector2(ttt[i], pm[j][i], alpha[j], pm[j][i]);
    }
}

static double
calc_pkp_mul ( fptype *restrict ttt
         , fptype *restrict cg_p
         , fptype msq_x4)
{
    fptype pkp = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:pkp)    
    for(int i = 0; i < num_floats_in_ks_array; ++i) {
        scalar_mult_add_su3_vector2(ttt[i], cg_p[i], msq_x4, ttt[i]); 
        pkp += cg_p[i] * ttt[i];
    }
    return double(pkp);
}

int
qphix_ks_multicg_offset ( void *t_src_arg
                        , void **t_dest_arg
                        , struct QuarkInvertControl *qic
                        , fptype *mass
                        , int num_offsets
                        , void *gll_arg[2]
                        , void *gfl_arg[2]
                        )
{
    //printf("qphix_ks_multicg_offset\n");
#if TIME_CG
    hrtimer_t t = hrtimer_get();
    CGPerfData cgPerfData;
#endif
    //printf("initializing ...\n");
    Gauge *gll[2];
    Gauge18 *gfl[2];
    KS *t_src       = (KS*)t_src_arg;
    KS *t_dest[MAXNDMASS];
    for(int i=0; i<num_offsets; ++i)
    t_dest[i]      = (KS*)t_dest_arg[i];
    gll[0] = (Gauge*) gll_arg[0];
    gll[1] = (Gauge*) gll_arg[1];

    gfl[0] = (Gauge18*) gfl_arg[0];

    int iteration;                      /* counter for iterations */
#ifdef FEWSUMS
    double actual_rsq;                  /* rsq from actual summation of resid */
    double c_tr,c_tt,tempsum[4];        /* Re<resid|ttt>, <ttt|ttt> */
#endif
    double rsq = 0, relrsq = 0;         /* resid**2, rel resid*2 */
    double oldrsq, pkp;                 /* last resid*2, pkp = cg_p.K.cg_p */
    double source_norm;                 /* squared magnitude of source vector */
    int nrestart;                       /* Restart counter */

    /* temporary storages */
    KS *ttt, *ttt_temp, /**cg_p,*/ *resid;
    int j, j_low;
    fptype * restrict shifts /*, shift0*/;
    fptype * restrict zeta_i, * restrict zeta_im1, * restrict zeta_ip1;
    fptype * restrict beta_i, * restrict beta_im1, * restrict alpha;
    int * restrict finished;
    KS ** restrict pm;
    int num_offsets_now;

    /* Unpack structure */
    int niter          = inv_arg->max;
    int max_restarts   = inv_arg->nrestart;
    fptype rsqmin      = resid_arg->resid;
    fptype relrsqmin   = resid_arg->relresid;
    int max_cg = max_restarts * niter;  /* Maximum number of iterations */
    int parity0 = inv_arg->parity, parity, otherparity; 
    /* if we want both parities, we will do even first. */
    switch(parity0){
    case(QPHIX_EVEN): parity=QPHIX_EVEN; otherparity=QPHIX_ODD; break;
    case(QPHIX_ODD):  parity=QPHIX_ODD; otherparity=QPHIX_EVEN; break;
    case(QPHIX_EVENODD):  parity=QPHIX_EVEN; otherparity=QPHIX_ODD; break;
    }

    /* Allocate the temps ttt, cg_p, resid */
    ttt      = (KS*)allocKS();
    ttt_temp = (KS*)allocKS();
    resid    = (KS*)allocKS();

    shifts = (fptype * )malloc(num_offsets*sizeof(fptype));
    zeta_i = (fptype * )malloc(num_offsets*sizeof(fptype));
    zeta_im1 = (fptype * )malloc(num_offsets*sizeof(fptype));
    zeta_ip1 = (fptype * )malloc(num_offsets*sizeof(fptype));
    beta_i = (fptype * )malloc(num_offsets*sizeof(fptype));
    beta_im1 = (fptype * )malloc(num_offsets*sizeof(fptype));
    alpha = (fptype * )malloc(num_offsets*sizeof(fptype));

    finished = (int *)malloc(sizeof(int)*num_offsets);

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

    pm = (KS **)malloc(num_offsets*sizeof(KS *));
    for(j=0; j<num_offsets; ++j)
        pm[j] = (KS *)allocKS();

  start:

    for(j = 0; j < num_offsets; j++)
    {
        finished[j] = 0;
        beta_im1[j] = -1.;
        zeta_im1[j] = 1.;
        zeta_i[j] = 1.;
        alpha[j] = 0.;
    }
    /* Source norm */
    source_norm = calc_source_norm((fptype*)t_src /*+ (parity == ODD ? num_floats_in_ks_array: 0)*/);
    g_doublesum( &source_norm );

    //printf("copy sources ...\n");
    /* copy sources to residues */
    for(j=0; j<num_offsets; ++j) if(j!=j_low)
        memcpy(pm[j], t_src, sizeof(fptype)*num_floats_in_ks_array);


    /* Start CG iterations */
    nrestart      = 0;
    iteration     = 0;
    num_offsets_now = num_offsets;
    resid_arg->size_r    = 0;
    resid_arg->size_relr = 0;

    //printf("start while loop...\n");
    while(true) {
        if( iteration % niter == 0 ||
            (rsqmin <= 0 || rsqmin > rsq/source_norm) &&
            (relrsqmin <= 0 || relrsqmin > relrsq) )
        {
#if TIME_CG
            cgPerfData.dslash_time -= hrtimer_get();
#endif
            //printf("calling qphix_ks_dslash ...\n");
            qphix_ks_dslash(t_dest[j_low]
                          , gll[otherparity]
                          , gfl[otherparity]
                          , ttt_temp
                          , otherparity);

            qphix_ks_dslash(ttt_temp
                          , gll[parity]
                          , gfl[parity]
                          , ttt
                          , parity);
#if TIME_CG
            cgPerfData.dslash_time += hrtimer_get();
#endif

#if TIME_CG
            cgPerfData.calc_rsq_time -= hrtimer_get();
#endif
            //printf("calling calc_rsq...\n");
            /* resid = t_src - (ttt + offset_low * t_dest) */
            /* pm = resid */
            rsq = calc_rsq_mul( (fptype*)ttt
                          , (fptype*)t_dest[j_low]
                          , (fptype*)t_src
                          , (fptype*)resid
                          , (fptype*)pm[j_low]
                          , offset_low
                          );
#if TIME_CG
            cgPerfData.calc_rsq_time += hrtimer_get();
#endif
#ifdef FEWSUMS
            actual_rsq = rsq; /* not yet summed over nodes */
#endif
            g_doublesum( &rsq );

            if(relrsqmin > 0) {
#if TIME_CG
                cgPerfData.relative_residue_time -= hrtimer_get();
#endif
                relrsq = relative_residue(resid, t_dest[j_low], parity);
#if TIME_CG
                cgPerfData.relative_residue_time += hrtimer_get();
#endif
            }

            resid_arg->final_rsq     = (fptype)rsq/source_norm;
            resid_arg->final_relrsq  = (fptype)relrsq;
            resid_arg->size_r        = resid_arg->final_rsq;
            resid_arg->size_relr     = resid_arg->final_relrsq;
            resid_arg->final_iters   = iteration;

            /* Converged */
            if( (rsqmin <= 0 || rsqmin > resid_arg->final_rsq) &&
                (relrsqmin <= 0 || relrsqmin > resid_arg->final_relrsq) ) {

                if(parity0==QPHIX_EVENODD && parity==QPHIX_EVEN) {
                    parity  = QPHIX_ODD;
                    otherparity = QPHIX_EVEN;
                    goto start;
                }
                else {
#if TIME_CG
                extern int myRank;
                if(myRank == 0)
                {
                        double cg_time = (hrtimer_get() - t)/1.0e09;
                        std::cout << myname << ": Total CG Time (s) : " << cg_time << " , Gflops : " << info->final_flop/cg_time << "\n";
                }
#if BW_CG
                cgPerfData.iters++;
                print_cg_timings(cgPerfData);
#endif
#endif
                     break;
                }
            }
            /* Otherwise */
            else {
                if( iteration >= max_cg || nrestart  >= max_restarts ) {
#if TIME_CG
                extern int myRank;
                if(myRank == 0)
                {
                        double cg_time = (hrtimer_get() - t)/1.0e09;
                        std::cout << myname << ": Total CG Time (s) : " << cg_time << " , Gflops : " << info->final_flop/cg_time << "\n";
                }
#if BW_CG
                print_cg_timings(cgPerfData);
#endif
#endif
                     break;
                }
            }

            nrestart++;

        } /* End checking residure */

        iteration++;

        oldrsq = rsq;

#if TIME_CG
        cgPerfData.dslash_time -= hrtimer_get();
#endif
        qphix_ks_dslash ( pm[j_low]
                        , gll[otherparity]
                        , gfl[otherparity]
                        , ttt_temp
                        , otherparity );

        qphix_ks_dslash ( ttt_temp
                        , gll[parity]
                        , gfl[parity]
                        , ttt
                        , parity );
#if TIME_CG
            cgPerfData.dslash_time += hrtimer_get();
#endif

#if TIME_CG
        cgPerfData.calc_pkp_time -= hrtimer_get();
#endif
        //printf("calling calc_pkp...\n");
        /* pkp = pm * (ttt + offset_low * pm) */
        /* ttt = ttt + offset_low * pm */
        pkp = calc_pkp_mul((fptype*)ttt, (fptype*)pm[j_low], offset_low);
#if TIME_CG
        cgPerfData.calc_pkp_time += hrtimer_get();
#endif
        g_doublesum( &pkp );

        /* Calculate params */
        zeta_ip1[j_low] = 1.0;
        beta_i[j_low] = -1.0 * rsq / pkp;
        for(j=0; j<num_offsets_now; ++j) if(j!=j_low) {
            fptype c = beta_i[j_low] * alpha[j_low] * (zeta_im1[j] - zeta_i[j]) +
                        zeta_im1[j] * beta_im1[j_low] * (1.0-shifts[j]*beta_i[j_low]);
            if(c==0.0) {
                zeta_ip1[j] = 0.0;
                finished[j] = 1;
            }
            else {
                zeta_ip1[j] = zeta_i[j] * zeta_im1[j] * beta_im1[j_low] / c;
                if(zeta_i[j] == 0.0)
                  finished[j] = 1;
                else
                  beta_i[j] = beta_i[j_low] * zeta_ip1[j] / zeta_i[j];
                if(beta_i[j]==0.0)
                  finished[j] = 1;
            }

            if( finished[j] == 1 && j == num_offsets_now-1 )
                num_offsets_now--;
          }

#if TIME_CG
          cgPerfData.calc_rsq2_time -= hrtimer_get();
#endif
          //printf("calling calc_rsq2_mul...\n");
          /* Update residue and solutions */
          /* t_dest = t_dest - beta_i * pm */
          /* resid = resid + beta_i * ttt */
          rsq = calc_rsq2_mul( (fptype **)t_dest
                         , (fptype *)beta_i
                         , (const fptype **)pm
                         , (fptype *)resid
                         , (fptype *)ttt
                         , num_offsets_now
                         , j_low
                         );
#if TIME_CG
          cgPerfData.calc_rsq2_time += hrtimer_get();
#endif

          if(relrsqmin > 0) {
#if TIME_CG
            cgPerfData.relative_residue_time -= hrtimer_get();
#endif
            relrsq = relative_residue(resid, t_dest[j_low], parity);
#if TIME_CG
            cgPerfData.relative_residue_time += hrtimer_get();
#endif
          }

          g_doublesum(&rsq);
/*
 *           resid_arg->size_r    = (fptype)rsq/source_norm;
 *                     resid_arg->size_relr = (fptype)relrsq;
 *                     */
          alpha[j_low] = (fptype)rsq / oldrsq;
          for(j=0; j<num_offsets_now; ++j) if(j!=j_low) {
                if(zeta_im1[j] != 0.0 && beta_im1[j_low] != 0.0)
                    alpha[j] = alpha[j_low] * zeta_i[j] * beta_im1[j] / (zeta_im1[j] * beta_im1[j_low]);
                else
                    alpha[j] = 0.0;

                if(alpha[j] == 0.0) {
                    finished[j] = 1;
                    if( j == num_offsets_now-1 )
                        num_offsets_now--;
                }
          }

          //printf("calling update_pm...\n");
          /* Update pm */
          /* pm = zeta_ip1 * resid + alpha * pm */
          update_pm( (fptype **)pm
                   , (fptype *)zeta_ip1
                   , (fptype *)resid
                   , (fptype *)alpha
                   , (fptype *)ttt
                   , num_offsets_now
                   );

          for(j=0; j<num_offsets_now; ++j)
          {
                zeta_im1[j] = zeta_i[j];
                zeta_i[j] = zeta_ip1[j];
                beta_im1[j] = beta_i[j];
          }

    } /* End while loop */

<<<<<<< HEAD
    if(myRank == 0) if(nrestart == max_restarts || iteration == max_cg) {
=======
    if(myRank==0) if(nrestart == max_restarts || iteration == max_cg) {
>>>>>>> gauge_force
        if(omp_get_thread_num() == 0) {
            std::cout << "qphix_congrad: CG not converged after " << iteration
                      << " iterations and " << nrestart << " restarts, \n"
                      << "rsq. = " << std::scientific << resid_arg->final_rsq
                      << " wanted " << rsqmin << " relrsq = "
                      << resid_arg->final_relrsq << " wanted " << relrsqmin
                      << std::endl;
        }

    }

    /* Free */
    freeKS(ttt);
    freeKS(ttt_temp);
    freeKS(resid);
    for(j=0; j<num_offsets; ++j)
        freeKS(pm[j]);
    free(pm);
    free(zeta_i);
    free(zeta_im1);
    free(zeta_ip1);
    free(beta_i);
    free(beta_im1);
    free(alpha);

#if TIME_CG
    extern int myRank;
    if(myRank == 0) std::cout << "Total CG Time (s) : " << (hrtimer_get() - t)/1.0e09 << "\n";
#endif
    return iteration;

}

#endif
