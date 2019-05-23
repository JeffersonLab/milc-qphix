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
  void *mpi_comm;                         /* Preinitialize MPI Communicator */
} QPHIX_layout_t;
#define QPHIX_LAYOUT_ZERO ((QPHIX_layout_t){NULL,NULL,0,NULL,0,NULL,0,0,0,NULL})

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

#define QPHIX_MAX_NAIK 20 /*verify this value */
typedef struct
{
  int n_naiks;
  double eps_naik[QPHIX_MAX_NAIK];
  //QPHIX_hisq_unitarize_group_t ugroup;
  //QPHIX_hisq_unitarize_method_t umethod;
  int ugroup;
  int umethod;
  double fat7_one_link;
  double fat7_three_staple;
  double fat7_five_staple;
  double fat7_seven_staple;
  double fat7_lepage;
  double asqtad_one_link;
  double asqtad_three_staple;
  double asqtad_five_staple;
  double asqtad_seven_staple;
  double asqtad_lepage;
  double asqtad_naik;
  double difference_one_link;
  double difference_naik;
} QPHIX_hisq_coeffs_t;

// Enumerations
typedef enum
{
  QPHIX_UNITARIZE_U3 = 0,
  QPHIX_UNITARIZE_SU3 = 1,
} QPHIX_hisq_unitarize_group_t;

typedef enum
{
  QPHIX_UNITARIZE_NONE = 0,
  QPHIX_UNITARIZE_RATIONAL = 3,
  QPHIX_UNITARIZE_ANALYTIC = 5,
} QPHIX_hisq_unitarize_method_t;

#define QPHIX_EPS_NAIK_ZERO {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
#define QPHIX_HISQ_COEFFS_ZERO \
  ((QPHIX_hisq_coeffs_t){1, QPHIX_EPS_NAIK_ZERO, QPHIX_UNITARIZE_U3, \
   QPHIX_UNITARIZE_RATIONAL, 0,0,0,0,0, 0,0,0,0,0,0, 0,0})

typedef struct
{
  char *tag;
  double value;
  void *extra;
} QPHIX_opt_t;

typedef struct
{
  int inited;
  int want_deps;
  int want_aux;
  QPHIX_hisq_unitarize_method_t umethod;
  int reunit_allow_svd;
  int reunit_svd_only;
  double reunit_svd_rel_error;
  double reunit_svd_abs_error;
  int svd_values_info;
  int use_fat7_lepage;
} QPHIX_hisq_links_t;

typedef struct
{
  double one_link;
  double three_staple;
  double five_staple;
  double seven_staple;
  double lepage;
  double naik;
} QPHIX_asqtad_coeffs_t;

typedef struct
{
  int inited;
  int fnmat_src_min;
  int veclength;
  double force_filter;
} QPHIX_hisq_force_t;

extern int qphix_even_sites_on_node;
extern QPHIX_hisq_links_t QPHIX_hisq_links;

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
QPHIX_status_t QPHIX_finalize();
QPHIX_status_t QPHIX_finalize_F();
QPHIX_status_t QPHIX_finalize_D();

#ifdef __cplusplus
}
#endif

#define QPHIX_OPP_PAR(p) ((p)==QPHIX_EVENODD ? (p) : (p)^11)
#define QPHIX_OPP_DIR(d) ((d)&1 ? (d)-1 : (d)+1)
#endif /* _QPHIX_INT_H */
