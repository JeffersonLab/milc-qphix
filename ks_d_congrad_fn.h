#ifndef _KS_D_CONGRAD_FN_H_
#define _KS_D_CONGRAD_FN_H_

#define MAXNDMASS 150

#ifdef __cplusplus
extern "C" {
#endif

#include "ks_config.h"   /* fptype */

/* Structure defining quark inversion parameters for most inverters */
struct QuarkInvertControl
{
    int prec;            /* Precision of the inversion 1 = single; 2 = double */
    int min;             /* Minimum number of iterations (being phased out) */
    int max;             /* Maximum number of iterations per restart */
    int nrestart;        /* Maximum restarts */
    int parity;          /* EVEN, ODD, or EVENANDODD (for some inverters) */
    int start_flag;      /* 0: use a zero initial guess; 1: use dest */
    int nsrc;            /* Number of source vectors */
    fptype resid;        /* Desired residual - normalized as :
                            sqrt(r*r)/sqrt(src_e*src_e) */
    fptype relresid;     /* Desired relative residual */
    fptype final_rsq;    /* Final true (absolute) residual */
    fptype final_relrsq; /* Final relative residual */
    fptype size_r;       /* Resulting cumulative residual */
    fptype size_relr;    /* Resulting cumulative relative residual */
    int converged;       /* Returned 0 if not converged; 1 if converged */
    int final_iters;     /* Number of iterations */
                         /* Add further parameters as needed...  */

};
    //qphix_quark_invert_control_t;

int 
qphix_ks_congrad_parity ( void *t_src_arg
                        , void *t_dest_arg
                        , struct QuarkInvertControl *qic
                        , fptype mass
                        , void *gll_arg[2]
                        , void *gfl_arg[2]
                        );

int
qphix_ks_congrad_two_src2 ( void *t_src1_arg
                          , void *t_src2_arg
                          , void *t_dest1_arg
                          , void *t_dest2_arg
                          , fptype mass1
                          , fptype mass2
                          , struct QuarkInvertControl *qic
                          , void *gll_arg[2]
                          , void *gfl_arg[2]
                          );

int
qphix_ks_multicg_offset ( void *t_src_arg
                        , void **t_dest_arg
                        , struct QuarkInvertControl *qic
                        , fptype *mass
                        , int num_offsets
                        , void *gll_arg[2]
                        , void *gfl_arg[2]
                        );

#ifdef __cplusplus
}
#endif

#endif
