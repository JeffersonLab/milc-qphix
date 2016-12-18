#ifndef _QPHIX_INT_H
#define _QPHIX_INT_H

/* This file defines miscellaneous types that may be used in all compilations */
/* It is included in qphix.h so it doesn't need to be included separately */

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  QPHIX_SUCCESS = 0,
  QPHIX_FAIL = 1,
  QPHIX_MEM_ERROR, /* memory errors such as failure to allocate, etc */
} QPHIX_status_t;

typedef enum {
  QPHIX_EVEN = 0x02,
  QPHIX_ODD = 0x01,
  QPHIX_EVENODD = 0x03
} QPHIX_evenodd_t;

typedef float QPHIX_F_Real;
typedef double QPHIX_D_Real;

typedef struct {
  int (*node_number)(const int coords[]); /* node no for given latt coord */
  int (*node_index)(const int coords[]);  /* site rank for given latt coord */
  int latdim;                             /* number of lattice dimensions */
  int *latsize;                           /* physical lattice lengths */
  int machdim;                            /* number of logical machine dims */
  int *machsize;                          /* logical grid lengths */
  int this_node;                          /* lexicographic node number */
  int sites_on_node;
  int even_sites_on_node;                 /* If needed */
} QPHIX_layout_t;
#define QPHIX_LAYOUT_ZERO ((QPHIX_layout_t){NULL,NULL,0,NULL,0,NULL,0,0})

typedef struct {
  double final_sec;        /* (out) total number of seconds used */
  double final_flop;       /* (out) total number of flops performed */
  QPHIX_status_t status;     /* (out) error status */
  int count1, count2;      /* (out) generic counters */
} QPHIX_info_t;
#define QPHIX_INFO_ZERO ((QPHIX_info_t){0,0,QPHIX_SUCCESS,0,0})

  /* these are quantities that apply to the gauge force calculation */
typedef struct {
   double plaquette;
   double rectangle;
   double parallelogram;
   double adjoint_plaquette;
} QPHIX_gauge_coeffs_t;
#define QPHIX_GAUGE_COEFFS_ZERO ((QPHIX_gauge_coeffs_t){0.,0.,0.,0.})

  /* these are quantities that apply to all masses in the multi inverter */
typedef struct {
  QPHIX_evenodd_t parity;
  int max;                  /* (in) max number of iterations */
  int restart;              /* (in) number of iterations before restart */
  int nrestart;             /* (in) number of restarts allowed */
  QPHIX_evenodd_t evenodd;  /* (in) subset of source vector */
} QPHIX_invert_arg_t;
#define QPHIX_INVERT_ARG_DEFAULT ((QPHIX_invert_arg_t){2000,1000,5,QPHIX_EVENODD})

  /* these are quantities that vary for each mass in the multi inverter */
typedef struct {
  double resid;          /* (in) desired squared residual. Ignored if 0 */
  double final_rsq;      /* (out) actual squared residual */
  double relresid;       /* (in) desired squared relative norm. Ignored if 0 */
  double final_rel;      /* (out) actual squared relative norm */
  double size_r;         /* resulting cumulative residual. Same
			    normalization as final_rsq. */
  double size_relr;      /* resulting cumulative relative
			    residual. Same normalization as
			    final_rel. */
  int final_iter;        /* (out) number of iterations done */
  int final_restart;     /* (out) number of restarts done */
} QPHIX_resid_arg_t;
#define QPHIX_RESID_ARG_DEFAULT ((QPHIX_resid_arg_t){1e-6,0,0,0,0,0})

  /**********************/
  /*  General routines  */
  /**********************/

QPHIX_status_t QPHIX_init_F(QPHIX_layout_t *layout);
QPHIX_status_t QPHIX_init_D(QPHIX_layout_t *layout);
QPHIX_status_t QPHIX_init(QPHIX_layout_t *layout);
QPHIX_status_t QPHIX_finalize(void);

#ifdef __cplusplus
}
#endif

#define QPHIX_OPP_PAR(p) ((p)==QPHIX_EVENODD ? (p) : (p)^11)
#define QPHIX_OPP_DIR(d) ((d)&1 ? (d)-1 : (d)+1)
#endif /* _QPHIX_INT_H */
